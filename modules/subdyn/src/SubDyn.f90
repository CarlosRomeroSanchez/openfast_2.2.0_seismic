!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.   
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
!> SubDyn is a time-domain structural-dynamics module for multi-member fixed-bottom substructures.
!! SubDyn relies on two main engineering schematizations: (1) a linear frame finite-element beam model (LFEB), and 
!! (2) a dynamics system reduction via Craig-Bampton�s (C-B) method, together with a Static-Improvement method, greatly reducing 
!!  the number of modes needed to obtain an accurate solution.   
Module SubDyn
   
   USE NWTC_Library
   USE NWTC_LAPACK
   USE SubDyn_Types
   USE SubDyn_Output
   USE SD_FEM
   
   IMPLICIT NONE

   PRIVATE
   
   !............................
   ! NOTE: for debugging, add preprocessor definition SD_SUMMARY_DEBUG
   !       this will add additional matrices to the SubDyn summary file.
   !............................
   TYPE(ProgDesc), PARAMETER  :: SD_ProgDesc = ProgDesc( 'SubDyn', '', '' )
      
   ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: SD_Init                           ! Initialization routine
   PUBLIC :: SD_End                            ! Ending routine (includes clean up)
   PUBLIC :: SD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
   PUBLIC :: SD_CalcOutput                     ! Routine for computing outputs
   PUBLIC :: SD_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   
CONTAINS

SUBROUTINE CreateTPMeshes( TP_RefPoint, inputMesh, outputMesh, ErrStat, ErrMsg )
   REAL(ReKi),                INTENT( IN    ) :: TP_RefPoint(3)
   TYPE(MeshType),            INTENT( INOUT ) :: inputMesh
   TYPE(MeshType),            INTENT( INOUT ) :: outputMesh
   INTEGER(IntKi),            INTENT(   OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   ! NOTE: The initialization of the fields for these meshes is to be handled by FAST/Driver
   CALL MeshCreate( BlankMesh        = inputMesh         &
                  ,IOS               = COMPONENT_INPUT   &
                  ,Nnodes            = 1                 &
                  ,ErrStat           = ErrStat           &
                  ,ErrMess           = ErrMsg            &
                  ,TranslationDisp   = .TRUE.            &
                  ,Orientation       = .TRUE.            &
                  ,TranslationVel    = .TRUE.            &
                  ,RotationVel       = .TRUE.            &
                  ,TranslationAcc    = .TRUE.            &
                  ,RotationAcc       = .TRUE.            )
   
   ! Create the node on the mesh
   CALL MeshPositionNode (   inputMesh           &
                           , 1                   &
                           , TP_RefPoint         &
                           , ErrStat             &
                           , ErrMsg              ) !note: assumes identiy matrix as reference orientation
   IF ( ErrStat >= AbortErrLev ) RETURN
      
   ! Create the mesh element
   CALL MeshConstructElement (   inputMesh          &
                               , ELEMENT_POINT      &                         
                               , ErrStat            &
                               , ErrMsg             &
                               , 1                  )
   CALL MeshCommit ( inputMesh, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   ! Create the Transition Piece reference point output mesh as a sibling copy of the input mesh
   CALL MeshCopy ( SrcMesh      = inputMesh              &
                  ,DestMesh     = outputMesh             &
                  ,CtrlCode     = MESH_SIBLING           &
                  ,IOS          = COMPONENT_OUTPUT       &
                  ,ErrStat      = ErrStat                &
                  ,ErrMess      = ErrMsg                 &
                  ,Force        = .TRUE.                 &
                  ,Moment       = .TRUE.                 ) 
END SUBROUTINE CreateTPMeshes

SUBROUTINE CreateY2Meshes( NNode, Nodes, NNodes_I, IDI, NNodes_L, IDL, NNodes_C, IDC, inputMesh, outputMesh, ErrStat, ErrMsg )
   INTEGER(IntKi),            INTENT( IN    ) :: NNode                     !total number of nodes in the structure, used to size the array Nodes, i.e. its rows
   REAL(ReKi),                INTENT( IN    ) :: Nodes(NNode, JointsCol)
   INTEGER(IntKi),            INTENT( IN    ) :: NNodes_I                  ! number interface nodes   i.e. Y2 stuff at the beginning
   INTEGER(IntKi),            INTENT( IN    ) :: IDI(NNodes_I*6)
   INTEGER(IntKi),            INTENT( IN    ) :: NNodes_L                  ! number interior nodes  (no constraints) i.e. Y2 stuff after interface stuff
   INTEGER(IntKi),            INTENT( IN    ) :: IDL(NNodes_L*6)
   INTEGER(IntKi),            INTENT( IN    ) :: NNodes_C                  ! number base reaction nodes  i.e. Y2 stuff after interior stuff
   INTEGER(IntKi),            INTENT( IN    ) :: IDC(NNodes_C*6)
   TYPE(MeshType),            INTENT( INOUT ) :: inputMesh
   TYPE(MeshType),            INTENT( INOUT ) :: outputMesh
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat                   ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg                    ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER         :: I                 ! generic counter variable
   INTEGER         :: nodeIndx
   
   CALL MeshCreate( BlankMesh        = inputMesh                           &
                  ,IOS               = COMPONENT_INPUT                     &
                  ,Nnodes            = NNodes_I + NNodes_L + NNodes_C      &
                  ,ErrStat           = ErrStat                             &
                  ,ErrMess           = ErrMsg                              &
                  ,Force             = .TRUE.                              &
                  ,Moment            = .TRUE.                              )
   !---------------------------------------------------------------------
   !    Interface nodes
   !---------------------------------------------------------------------
   DO I = 1,NNodes_I 
      ! Create the node on the mesh
      nodeIndx = IDI(I*6) / 6     !integer division gives me the actual node index, is it true? Yes it is not the nodeID
      CALL MeshPositionNode (   inputMesh           &
                              , I                   &
                              , Nodes(nodeIndx,2:4) &  ! position
                              , ErrStat             &
                              , ErrMsg              )
      IF ( ErrStat /= ErrID_None ) RETURN

      ! Create the mesh element
      CALL MeshConstructElement (   inputMesh          &
                                  , ELEMENT_POINT      &                         
                                  , ErrStat            &
                                  , ErrMsg             &
                                  , I                  )
   END DO
   
   !---------------------------------------------------------------------
   !    Interior nodes
   !---------------------------------------------------------------------
   DO I = 1,NNodes_L 
      ! Create the node on the mesh
      nodeIndx = IDL(I*6) / 6     !integer division gives me the actual node index, is it true? Yes it is not the nodeID of the input file that may not be sequential, but the renumbered list of nodes
      CALL MeshPositionNode (   inputMesh           &
                              , I + NNodes_I        &
                              , Nodes(nodeIndx,2:4) &
                              , ErrStat             &
                              , ErrMsg              )
      IF ( ErrStat /= ErrID_None ) RETURN

      ! Create the mesh element
      CALL MeshConstructElement (   inputMesh          &
                                  , ELEMENT_POINT      &                         
                                  , ErrStat            &
                                  , ErrMsg             &
                                  , I + NNodes_I       )
   END DO
   
   !---------------------------------------------------------------------
   !    Base Reaction nodes
   !---------------------------------------------------------------------
   DO I = 1,NNodes_C 
      ! Create the node on the mesh
      nodeIndx = IDC(I*6) / 6     !integer division gives me the actual node index, is it true? Yes it is not the nodeID
      CALL MeshPositionNode (   inputMesh                 &
                              , I + NNodes_I + NNodes_L   &
                              , Nodes(nodeIndx,2:4)       &  
                              , ErrStat                   &
                              , ErrMsg                    )
      IF ( ErrStat /= ErrID_None ) RETURN
      
      ! Create the mesh element
      CALL MeshConstructElement (   inputMesh                 &
                                  , ELEMENT_POINT             &                         
                                  , ErrStat                   &
                                  , ErrMsg                    &
                                  , I + NNodes_I + NNodes_L   )
   END DO
   CALL MeshCommit ( inputMesh, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
         
   ! Create the Interior Points output mesh as a sibling copy of the input mesh
   CALL MeshCopy (    SrcMesh      = inputMesh              &
                     ,DestMesh     = outputMesh             &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat                &
                     ,ErrMess      = ErrMsg                 &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.            ) 
   
    ! Set the Orientation (rotational) field for the nodes based on assumed 0 (rotational) deflections
    !Identity should mean no rotation, which is our first guess at the output -RRD
    CALL Eye( outputMesh%Orientation, ErrStat, ErrMsg )         
        
END SUBROUTINE CreateY2Meshes
!------------------------------------------------------------------------------------------------------
!> Set the index array that maps SD internal nodes to the Y2Mesh nodes.
!! NOTE: SDtoMesh is not checked for size, nor are the index array values checked for validity, 
!!       so this routine could easily have segmentation faults if any errors exist.
SUBROUTINE SD_Y2Mesh_Mapping(p, SDtoMesh )
   TYPE(SD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   INTEGER(IntKi),               INTENT(  OUT)  :: SDtoMesh(:) !< index/mapping of mesh nodes with SD mesh
   ! locals
   INTEGER(IntKi)                               :: i
   INTEGER(IntKi)                               :: SDnode
   INTEGER(IntKi)                               :: y2Node

   y2Node = 0
   ! Interface nodes (IDI)
   DO I = 1,SIZE(p%IDI,1)/6
      y2Node = y2Node + 1      
      SDnode = p%IDI(I*6) / 6     !integer division gives me the actual node index; it is not the nodeID
      SDtoMesh( SDnode ) = y2Node ! TODO add safety check
   END DO
   
   ! Interior nodes (IDL)
   DO I = 1,SIZE(p%IDL,1)/6 
      y2Node = y2Node + 1      
      SDnode = p%IDL(I*6) / 6     !integer division gives me the actual node index; it is not the nodeID
      SDtoMesh( SDnode ) = y2Node ! TODO add safety check
   END DO

   ! Base Reaction nodes (IDC)
   DO I = 1,SIZE(p%IDC,1)/6 
      y2Node = y2Node + 1      
      SDnode = p%IDC(I*6) / 6     !integer division gives me the actual node index; it is not the nodeID
      SDtoMesh( SDnode ) = y2Node ! TODO add safety check
   END DO

END SUBROUTINE SD_Y2Mesh_Mapping


!---------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE SD_Init( InitInput, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
   TYPE(SD_InitInputType),       INTENT(IN   )  :: InitInput   !< Input data for initialization routine         
   TYPE(SD_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
   TYPE(SD_ParameterType),       INTENT(  OUT)  :: p           !< Parameters
   TYPE(SD_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
   TYPE(SD_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
   TYPE(SD_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
   TYPE(SD_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states
   TYPE(SD_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                               !!    only the output mesh is initialized)
   REAL(DbKi),                   INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that
                                                               !!   (1) Mod1_UpdateStates() is called in loose coupling &
                                                               !!   (2) Mod1_UpdateDiscState() is called in tight coupling.
                                                               !!   Input is the suggested time from the glue code;
                                                               !!   Output is the actual coupling interval that will be used
                                                               !!   by the glue code.
   TYPE(SD_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables
   TYPE(SD_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! local variables
   TYPE(SD_InitType)    :: Init
   TYPE(CB_MatArrays)   :: CBparams      ! CB parameters to be stored and written to summary file
   TYPE(FEM_MatArrays)  :: FEMparams     ! FEM parameters to be stored and written to summary file
   INTEGER(IntKi)       :: ErrStat2      ! Error status of the operation
   CHARACTER(ErrMsgLen) :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   
   ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! Initialize the NWTC Subroutine Library
   CALL NWTC_Init( )

   ! Display the module information
   CALL DispNVD( SD_ProgDesc )   
   InitOut%Ver = SD_ProgDesc
   
   ! transfer glue-code information to data structure for SubDyn initialization:
   Init%g           = InitInput%g   
   Init%TP_RefPoint = InitInput%TP_RefPoint
   Init%SubRotateZ  = InitInput%SubRotateZ
   p%NAvgEls        = 2

   !bjj added this ugly check (mostly for checking SubDyn driver). not sure if anyone would want to play with different values of gravity so I don't return an error.
   IF (Init%g < 0.0_ReKi ) CALL ProgWarn( ' SubDyn calculations use gravity assuming it is input as a positive number; the input value is negative.' ) 
   
   ! Establish the GLUECODE requested/suggested time step.  This may be overridden by SubDyn based on the SDdeltaT parameter of the SubDyn input file.
   Init%DT  = Interval
   IF ( LEN_TRIM(Init%RootName) == 0 ) THEN
      CALL GetRoot( InitInput%SDInputFile, Init%RootName )
   ELSE
      Init%RootName = TRIM(InitInput%RootName)//'.SD'
   END IF
   
   ! Parse the SubDyn inputs 
   CALL SD_Input(InitInput%SDInputFile, Init, p, ErrStat2, ErrMsg2); if(Failed()) return
   
   ! Discretize the structure according to the division size 
   ! sets Init%NNode, Init%NElm
   CALL SD_Discrt(Init,p, ErrStat2, ErrMsg2); if(Failed()) return
      
   ! Assemble Stiffness and mass matrix
   CALL AssembleKM(Init,p, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Calculate values for FEMparams (for summary file output only
   ! Solve dynamics problem
   FEMparams%NOmega = Init%TDOF - p%Nreact*6 !removed an extra "-6"  !Note if fixity changes at the reaction points, this will need to change
     
   CALL AllocAry(FEMparams%Omega,            FEMparams%NOmega, 'FEMparams%Omega', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL AllocAry(FEMparams%Modes, Init%TDOF, FEMparams%NOmega, 'FEMparams%Modes', ErrStat2, ErrMsg2 ); if(Failed()) return
   
   ! We call the EigenSolver here only so that we get a print-out the eigenvalues from the full system (minus Reaction DOF)
   ! The results, Phi is not used in the remainder of this Init subroutine, Omega goes to outsummary.
   CALL EigenSolve( Init%K, Init%M, Init%TDOF, FEMparams%NOmega, .True., Init, p, FEMparams%Modes, FEMparams%Omega, ErrStat2, ErrMsg2 ); if(Failed()) return

   IF(Init%CBMod) THEN ! C-B reduction  

      ! --- Craig-Bampton reduction (sets many parameters)
      CALL Craig_Bampton(Init, p, CBparams, ErrStat2, ErrMsg2); if(Failed()) return

      ! --- Initial system states 
      IF ( p%qmL > 0 ) THEN
         CALL AllocAry(x%qm,       p%qmL, 'x%qm',       ErrStat2, ErrMsg2 ); if(Failed()) return
         CALL AllocAry(x%qmdot,    p%qmL, 'x%qmdot',    ErrStat2, ErrMsg2 ); if(Failed()) return
         CALL AllocAry(m%qmdotdot, p%qmL, 'm%qmdotdot', ErrStat2, ErrMsg2 ); if(Failed()) return
         x%qm      = 0.0_ReKi   
         x%qmdot   = 0.0_ReKi
         m%qmdotdot= 0.0_ReKi
      END IF

   ELSE  !Full FEM directly

      CALL FULLFEM(Init, p, FEMparams, Errstat2, ErrMsg2); if(Failed()) return

      ! --- Even if the state variables are now UL instead of qm, the name of the variable is kept the same
      !     in order to avoid changing a number of routines just because of the name of the variable
      CALL AllocAry(x%qm,       p%DOFL, 'x%qm',       ErrStat2, ErrMsg2 ); if(Failed()) return
      CALL AllocAry(x%qmdot,    p%DOFL, 'x%qmdot',    ErrStat2, ErrMsg2 ); if(Failed()) return
      CALL AllocAry(m%qmdotdot, p%DOFL, 'm%qmdotdot', ErrStat2, ErrMsg2 ); if(Failed()) return
      x%qm      = 0.0_ReKi   
      x%qmdot   = 0.0_ReKi
      m%qmdotdot= 0.0_ReKi

   ENDIF
   
   xd%DummyDiscState  = 0.0_ReKi
   z%DummyConstrState = 0.0_ReKi

   ! Allocate OtherState%xdot if using multi-step method; initialize n
   IF ( ( p%IntMethod .eq. 2) .OR. ( p%IntMethod .eq. 3)) THEN
      !bjj: note that the way SD_UpdateStates is implemented, "n" doesn't need to be initialized here
      Allocate( OtherState%xdot(4), STAT=ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating OtherState%xdot', ErrStat, ErrMsg, 'SD_Init' )
         CALL CleanUp()
         RETURN
      END IF
   ENDIF
 
   ! Allocate miscellaneous variables, used only to avoid temporary copies of variables allocated/deallocated and sometimes recomputed each time
   CALL AllocMiscVars(p, m, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Write the summary file
   IF ( Init%SSSum ) THEN 
      ! note p%KBB/MBB are KBBt/MBBt
      ! Write a summary of the SubDyn Initialization                     
      CALL OutSummary(Init,p,FEMparams,CBparams,  ErrStat2, ErrMsg2); if(Failed()) return
      IF( ALLOCATED(Init%K) ) DEALLOCATE(Init%K)
      IF( ALLOCATED(Init%C) ) DEALLOCATE(Init%C)
      IF( ALLOCATED(Init%M) ) DEALLOCATE(Init%M)     
   ENDIF 
      
   ! --- Initialize Inputs and Outputs
   ! Create the input and output meshes associated with Transition Piece reference point       
   CALL CreateTPMeshes( InitInput%TP_RefPoint, u%TPMesh, y%Y1Mesh, ErrStat2, ErrMsg2 ); if(Failed()) return
   
   ! Construct the input mesh for the interior nodes which result from the Craig-Bampton reduction
   CALL CreateY2Meshes( Init%NNode, Init%Nodes, Init%NInterf, p%IDI, p%NNodes_L, p%IDL, p%NReact, p%IDC, u%LMesh, y%Y2Mesh, ErrStat2, ErrMsg2 ); if(Failed()) return
   
   ! Initialize the outputs & Store mapping between nodes and elements  
   CALL SDOUT_Init( Init, y, p, m, InitOut, InitInput%WtrDpth, ErrStat2, ErrMsg2 ); if(Failed()) return
   
   ! Determine if we need to perform output file handling
   IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN  
       CALL SDOUT_OpenOutput( SD_ProgDesc, Init%RootName, p, InitOut, ErrStat2, ErrMsg2 ); if(Failed()) return
   END IF
      
   
   ! Tell GLUECODE the SubDyn timestep interval 
   Interval = p%SDdeltaT
   CALL CleanUp()

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Init') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed
   
   SUBROUTINE CleanUp()   
      CALL SD_DestroyInitType(Init,   ErrStat2, ErrMsg2)
      CALL SD_DestroyCB_MatArrays(  CBparams,  ErrStat2, ErrMsg2 )  ! local variables
      CALL SD_DestroyFEM_MatArrays( FEMparams, ErrStat2, ErrMsg2 )  ! local variables
   END SUBROUTINE CleanUp

END SUBROUTINE SD_Init

!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, discrete, constraint, and other states are updated for t + Interval.
SUBROUTINE SD_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      REAL(DbKi),                         INTENT(IN   ) :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
      TYPE(SD_InputType),                 INTENT(INOUT) :: Inputs(:)       !< Inputs at Times
      REAL(DbKi),                         INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(SD_ParameterType),             INTENT(IN   ) :: p               !< Parameters
      TYPE(SD_ContinuousStateType),       INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                           !!   Output: Continuous states at t + Interval
      TYPE(SD_DiscreteStateType),         INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                           !!   Output: Discrete states at t + Interval
      TYPE(SD_ConstraintStateType),       INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                           !!   Output: Constraint states at t + Interval
      TYPE(SD_OtherStateType),            INTENT(INOUT) :: OtherState      !< Input: Other states at t;
                                                                           !!   Output: Other states at t + Interval
      TYPE(SD_MiscVarType),               INTENT(INOUT) :: m               !< Misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None
      ! Initialize variables
      ErrStat   = ErrID_None           ! no error has occurred
      ErrMsg    = ""

      IF ( p%SeismicInp ) CALL InterpSeismicSignal(t, p, m)
      
      IF ((p%CBMod).AND.( p%qml == 0)) RETURN ! no retained modes = no states
        
      IF (p%IntMethod .eq. 1) THEN
         CALL SD_RK4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      ELSEIF (p%IntMethod .eq. 2) THEN
         CALL SD_AB4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      ELSEIF (p%IntMethod .eq. 3) THEN
         CALL SD_ABM4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      ELSE  
         CALL SD_AM2( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      END IF
      
END SUBROUTINE SD_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE SD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(SD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
      TYPE(SD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
      TYPE(SD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(SD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(SD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                  !!   nectivity information does not have to be recalculated)
      TYPE(SD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      !locals
      INTEGER(IntKi)               :: L1,L2       ! partial Lengths of state and input arrays
      INTEGER(IntKi)               :: I,J         ! Counters
      REAL(ReKi)                   :: AllOuts(0:MaxOutPts+p%OutAllInt*p%OutAllDims)
      REAL(ReKi)                   :: rotations(3)
      REAL(ReKi)                   :: ULS(p%DOFL),  UL0m(p%DOFL),  FLt(p%DOFL)  ! Temporary values in static improvement method
      REAL(ReKi)                   :: Y1(6)
      INTEGER(IntKi)               :: startDOF
      REAL(ReKi)                   :: DCM(3,3),junk(6,p%NNodes_L)
      REAL(ReKi)                   :: HydroForces(6*p%NNodes_I) !  !Forces from all interface nodes listed in one big array  ( those translated to TP ref point HydroTP(6) are implicitly calculated in the equations)
      TYPE(SD_ContinuousStateType) :: dxdt        ! Continuous state derivatives at t- for qmdotdot purposes only
      INTEGER(IntKi)               :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(ErrMsgLen)         :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None
      REAL(ReKi), ALLOCATABLE      :: seismic_force_temp(:)
      CHARACTER(*), PARAMETER      :: RoutineName = 'SD_CalcOutput'
                                                 
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

                                   
      ! Compute the small rotation angles given the input direction cosine matrix
      rotations  = GetSmllRotAngs(u%TPMesh%Orientation(:,:,1), ErrStat2, Errmsg2); if(Failed()) return
      
      ! Inputs at the transition piece:
      m%u_TP       = (/REAL(u%TPMesh%TranslationDisp(:,1),ReKi), rotations/)
      m%udot_TP    = (/u%TPMesh%TranslationVel( :,1), u%TPMesh%RotationVel(:,1)/)
      m%udotdot_TP = (/u%TPMesh%TranslationAcc( :,1), u%TPMesh%RotationAcc(:,1)/)
      ! Inputs on interior nodes:
      CALL ConstructUFL( u, p, m%UFL )

      !________________________________________
      ! Set motion outputs on y%Y2mesh
      !________________________________________
      ! Y2 = C2*x + D2*u + F2 

      m%UR_bar        =                                      matmul( p%TI      , m%u_TP       )  ! UR_bar         [ Y2(1) =       0*x(1) + D2(1,1)*u(1) ]      
      m%UR_bar_dot    =                                      matmul( p%TI      , m%udot_TP    )  ! UR_bar_dot     [ Y2(3) =       0*x(1) + D2(3,2)*u(2) ]
      m%UR_bar_dotdot =                                      matmul( p%TI      , m%udotdot_TP )  ! U_R_bar_dotdot [ Y2(5) =       0*x(2) + D2(5,3)*u(3) ] 

      IF(p%CBMod) THEN ! C-B reduction    

        IF ( p%qml > 0) THEN
         m%UL            = matmul( p%PhiM,  x%qm    )      + matmul( p%PhiRb_TI, m%u_TP       )  ! UL             [ Y2(2) = C2(2,1)*x(1) + D2(2,1)*u(1) ] : IT MAY BE MODIFIED LATER IF STATIC IMPROVEMENT
         m%UL_dot        = matmul( p%PhiM,  x%qmdot )      + matmul( p%PhiRb_TI, m%udot_TP    )  ! UL_dot         [ Y2(4) = C2(2,2)*x(2) + D2(4,2)*u(2) ]      
         
         m%UL_dotdot     = matmul( p%C2_61, x%qm    )   + matmul( p%C2_62   , x%qmdot )       &  ! UL_dotdot      [ Y2(6) = C2(6,1)*x(1) + C2(6,2)*x(2) ...
                         + matmul( p%D2_62, m%udot_TP )                                       &  !                        + D2(6,2)*u(2)
                         + matmul( p%D2_63, m%udotdot_TP ) + matmul( p%D2_64,    m%UFL      ) &  !                        + D2(6,3)*u(3) + D2(6,4)*u(4) ...  ! -> bjj: this line takes up a lot of time. are any matrices sparse?
                                  + p%F2_61

         IF (p%SeismicInp) THEN
            m%UL         =  m%UL        + matmul ( p%PhiRbase , p%RRbase ) * m%Ug 
            m%UL_dot     =  m%UL_dot    + matmul ( p%PhiRbase , p%RRbase ) * m%Udotg 
            m%UL_dotdot  =  m%UL_dotdot + matmul ( p%PhiM ,  m%Ug * p%FMIMKP + m%Udotg * p%FMIMCP - m%Udotg * p%FMIMMP + m%PHIg * p%FMIMKP_PHI + m%PHIdotg * p%FMIMCP_PHI - m%PHIdotg * p%FMIMMP_PHI + m%Vg * p%FMIMKP_V + m%Vdotg * p%FMIMCP_V - m%Vdotg * p%FMIMMP_V ) &
                                        + matmul ( p%PhiRbase , p%RRbase ) * m%Uddotg  
         ENDIF
                                                                                 !                        + F2(6) ]                  
        ELSE ! There are no states when p%qml=0 (i.e., no retained modes: p%Nmodes=0), so we omit those portions of the equations
         m%UL            =                                   matmul( p%PhiRb_TI, m%u_TP       )  ! UL             [ Y2(2) =       0*x(1) + D2(2,1)*u(1) ] : IT MAY BE MODIFIED LATER IF STATIC IMPROVEMENT
         m%UL_dot        =                                   matmul( p%PhiRb_TI, m%udot_TP    )  ! UL_dot         [ Y2(4) =       0*x(2) + D2(4,2)*u(2) ]      
         m%UL_dotdot     =                                   matmul( p%PhiRb_TI, m%udotdot_TP )  ! UL_dotdot      [ Y2(6) =       0*x(:) + D2(6,3)*u(3) + 0*u(4) + 0]

         IF (p%SeismicInp) THEN
            m%UL         =  m%UL        + matmul ( p%PhiRbase , p%RRbase ) * m%Ug 
            m%UL_dot     =  m%UL_dot    + matmul ( p%PhiRbase , p%RRbase ) * m%Udotg 
            m%UL_dotdot  =  m%UL_dotdot + matmul ( p%PhiRbase , p%RRbase ) * m%Uddotg  
         ENDIF

        END IF
      
      !STATIC IMPROVEMENT METHOD  ( modify UL )
        IF (p%SttcSolve) THEN
         FLt  = MATMUL(p%PhiL_T,                  m%UFL + p%FGL)  ! -> bjj: todo: this line takes up A LOT of time. is PhiL sparse???? no (solution: don't call this routine thousands of time to calculate the jacobian)
         ULS  = MATMUL(p%PhiLInvOmgL2,            FLt          )  ! -> bjj: todo: this line takes up A LOT of time. is PhiL sparse????
         m%UL = m%UL + ULS 
          
         IF ( p%qml > 0) THEN
            UL0M = MATMUL(p%PhiLInvOmgL2(:,1:p%qmL), FLt(1:p%qmL)       )
            m%UL = m%UL - UL0M 
         END IF          
        ENDIF   
      
      ELSE !Full FEM

        ! note that this re-sets m%udotdot_TP, m%dot_TP, m%u_TP and m%UFL
        !but they are the same values as earlier in this routine so it doesn't change results in SDOut_MapOutputs()
        CALL SD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat2, ErrMsg2 ); if(Failed()) return
        !Save the acceleration 
        m%qmdotdot=dxdt%qmdot

        m%UL        = x%qm 
        m%UL_dot    = x%qmdot
        m%UL_dotdot = m%qmdotdot  

      ENDIF
                                                      
      ! --------------------------------------------------------------------------------- 
      ! Place the outputs onto interface node portion of Y2 output mesh        
      ! ---------------------------------------------------------------------------------
      DO I = 1, p%NNodes_I 
         startDOF = (I-1)*6 + 1
         ! Construct the direction cosine matrix given the output angles
         CALL SmllRotTrans( 'UR_bar input angles', m%UR_bar(startDOF + 3), m%UR_bar(startDOF + 4), m%UR_bar(startDOF + 5), DCM, '', ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput')

         y%Y2mesh%TranslationDisp (:,I)     = m%UR_bar  (       startDOF     : startDOF + 2 )
         y%Y2mesh%Orientation     (:,:,I)   = DCM
         y%Y2mesh%TranslationVel  (:,I)     = m%UR_bar_dot (    startDOF     : startDOF + 2 )
         y%Y2mesh%RotationVel     (:,I)     = m%UR_bar_dot (    startDOF + 3 : startDOF + 5 )
         y%Y2mesh%TranslationAcc  (:,I)     = m%UR_bar_dotdot ( startDOF     : startDOF + 2 )
         y%Y2mesh%RotationAcc     (:,I)     = m%UR_bar_dotdot ( startDOF + 3 : startDOF + 5 )
                  
      ENDDO
     
      ! --------------------------------------------------------------------------------- 
      ! Place the outputs onto interior node portion of Y2 output mesh 
      ! ---------------------------------------------------------------------------------      
      DO I = 1, p%NNodes_L   !Only interior nodes here     
         ! starting index in the master arrays for the current node    
         startDOF = (I-1)*6 + 1
         
         ! index into the Y2Mesh
         J = p%NNodes_I + I
       
         ! Construct the direction cosine matrix given the output angles
         CALL SmllRotTrans( 'UL input angles', m%UL(startDOF + 3), m%UL(startDOF + 4), m%UL(startDOF + 5), DCM, '', ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput')
         
         ! Y2 = Interior node displacements and velocities  for use as inputs to HydroDyn
         y%Y2mesh%TranslationDisp (:,J)     = m%UL     ( startDOF     : startDOF + 2 )
         y%Y2mesh%Orientation     (:,:,J)   = DCM
         y%Y2mesh%TranslationVel  (:,J)     = m%UL_dot ( startDOF     : startDOF + 2 )
         y%Y2mesh%RotationVel     (:,J)     = m%UL_dot ( startDOF + 3 : startDOF + 5 )
         
      END DO
      
      !Repeat for the acceleration, there should be a way to combine into 1 loop
      L1 = p%NNodes_I+1
      L2 = p%NNodes_I+p%NNodes_L
      junk=   RESHAPE(m%UL_dotdot,(/6  ,p%NNodes_L/)) 
      y%Y2mesh%TranslationAcc (  :,L1:L2)   = junk(1:3,:) 
      y%Y2mesh%RotationAcc    (  :,L1:L2)   = junk(4:6,:) 
      
      ! ---------------------------------------------------------------------------------
      ! Base reaction nodes
      ! ---------------------------------------------------------------------------------
      L1 = p%NNodes_I+p%NNodes_L+1   
      L2 = p%NNodes_I+p%NNodes_L+p%NReact

      y%Y2mesh%TranslationDisp(  :,L1:L2)   = 0.0
      CALL Eye( y%Y2mesh%Orientation(:,:,L1:L2), ErrStat2, ErrMsg2 ) ; if(Failed()) return

      y%Y2mesh%TranslationVel (  :,L1:L2)   = 0.0
      y%Y2mesh%RotationVel    (  :,L1:L2)   = 0.0
      y%Y2mesh%TranslationAcc (  :,L1:L2)   = 0.0
      y%Y2mesh%RotationAcc    (  :,L1:L2)   = 0.0

      !________________________________________
      ! Set loads outputs on y%Y1Mesh
      !________________________________________
      ! ---------------------------------------------------------------------------------
      !Y1= TP reaction Forces, i.e. force that the jacket exerts onto the TP and above  
      ! ---------------------------------------------------------------------------------
      ! Eq. 15: Y1 = -(C1*x + D1*u + FY)  [note the negative sign!!!!]
      !NEED TO ADD HYDRODYNAMIC FORCES AT THE Interface NODES
        !Aggregate the forces and moments at the interface nodes to the reference point
        !TODO: where are these HydroTP, HydroForces documented?
      DO I = 1, p%NNodes_I 
         startDOF = (I-1)*6 + 1
         !Take care of Hydrodynamic Forces that will go into INterface Forces later
         HydroForces(startDOF:startDOF+5 ) =  (/u%LMesh%Force(:,I),u%LMesh%Moment(:,I)/)  !(6,NNODES_I)
      ENDDO

      IF(p%CBMod) THEN ! C-B reduction        
                
        !HydroTP =  matmul(transpose(p%TI),HydroForces) ! (6,1) calculated below
        ! note: matmul( HydroForces, p%TI ) = matmul( transpose(p%TI), HydroForces) because HydroForces is 1-D            
        IF ( p%qml > 0) THEN
         Y1 = -(   matmul(p%C1_11, x%qm) + matmul(p%C1_12,x%qmdot)                                    &  ! -(   C1(1,1)*x(1) + C1(1,2)*x(2)
                 + matmul(p%D1_12, m%udot_TP)                                                         &  !    + D1(1,2)*u(2) +
                 + matmul(p%KBB,   m%u_TP) + matmul(p%D1_13, m%udotdot_TP) + matmul(p%D1_14, m%UFL)   &  !    + D1(1,1)*u(1) + 0*u(2) + D1(1,3)*u(3) + D1(1,4)*u(4)
                 - matmul( HydroForces, p%TI )  + p%FY )                                                                            !    + D1(1,5)*u(5) + Fy(1) )

         IF (p%SeismicInp) THEN
           Y1 = Y1 - ( matmul(p%MBM ,  m%Ug * p%FMIMKP + m%Udotg * p%FMIMCP - m%Uddotg * p%FMIMMP + m%PHIg * p%FMIMKP_PHI + m%PHIdotg * p%FMIMCP_PHI - m%PHIddotg * p%FMIMMP_PHI  + m%Vg * p%FMIMKP_V + m%Vdotg * p%FMIMCP_V - m%Vddotg * p%FMIMMP_V)      &
                     - (matmul( TRANSPOSE(p%TI),m%Ug * p%FRIMKP + m%Udotg * p%FRIMCP - m%Uddotg * p%FRIMMP + m%PHIg * p%FRIMKP_PHI + m%PHIdotg * p%FRIMCP_PHI - m%PHIddotg * p%FRIMMP_PHI + m%Vg * p%FRIMKP_V + m%Vdotg * p%FRIMCP_V - m%Vddotg * p%FRIMMP_V )))
         ENDIF

        ELSE ! No retained modes, so there are no states
         Y1 = -( matmul(p%KBB,   m%u_TP) + matmul(p%CBB, m%udot_TP) + matmul(p%MBB, m%udotdot_TP)    &  ! -(  0*x + D1(1,1)*u(1) + D1(1,2)*u(2) + 0*u(2) + D1(1,3)*u(3) 
                + matmul(p%D1_14, m%UFL) - matmul( HydroForces, p%TI )  + p%FY )                        !   + D1(1,4)*u(4) + D1(1,5)*u(5) + Fy(1) )

         IF (p%SeismicInp) THEN
           Y1 = Y1 - ( - matmul( TRANSPOSE(p%TI) ,  m%Ug * p%FRIMKP + m%Udotg * p%FRIMCP - m%Uddotg * p%FRIMMP + m%PHIg * p%FRIMKP_PHI + m%PHIdotg * p%FRIMCP_PHI - m%PHIddotg * p%FRIMMP_PHI + m%Vg * p%FRIMKP_V + m%Vdotg * p%FRIMCP_V - m%Vddotg * p%FRIMMP_V ))
         ENDIF

        END IF

       ELSE !Full FEM

         Y1 = - (   matmul(p%C1_11 , x%qm)            &
                  + matmul(p%C1_12 , x%qmdot)         &
                  + matmul(p%D1_11 , m%u_TP)          &
                  + matmul(p%D1_12 , m%udot_TP)       &
                  + matmul(p%D1_13 , m%udotdot_TP)    &
                  + matmul(p%D1_14 , m%UFL)           &
                  - matmul( HydroForces, p%TI ) + p%FY )

         IF (p%SeismicInp) THEN
             CALL AllocAry(seismic_force_temp, p%DOFL, 'seismic_force_temp', ErrStat2, ErrMsg2 ); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput' )

             seismic_force_temp = m%Ug * matmul(p%D1_14 , p%FSISKL)
             Y1 = Y1 - seismic_force_temp

             seismic_force_temp = m%Udotg * matmul(p%D1_14 , p%FSISCL)
             Y1 = Y1 - seismic_force_temp

             seismic_force_temp = m%Uddotg * matmul(p%D1_14 , p%FSISML)
             Y1 = Y1 + seismic_force_temp

             Y1 = Y1 + p%FSISKI*m%Ug + p%FSISCI*m%Udotg + p%FSISMI*m%Uddotg

             DEALLOCATE(seismic_force_temp)
         ENDIF

       ENDIF
      
      ! values on the interface mesh are Y1 (SubDyn forces) + Hydrodynamic forces
      y%Y1Mesh%Force (:,1) = Y1(1:3) 
      y%Y1Mesh%Moment(:,1) = Y1(4:6) 
            
     !________________________________________
     ! CALCULATE OUTPUT TO BE WRITTEN TO FILE 
     !________________________________________
     ! OutSwtch determines whether or not to actually output results via the WriteOutput array
     !    0 = No one needs the SubDyn outputs provided via the WriteOutput array.
     !    1 = SubDyn will generate an output file of its own.  
     !    2 = the caller will handle the outputs, but SubDyn needs to provide them.
     !    3 = Both 1 and 2
      IF ( p%OutSwtch > 0 ) THEN
         ! call CalcContStateDeriv one more time to store these qmdotdot for debugging purposes in the output file
         !find xdot at t
         IF ( p%NModes > 0 ) THEN
            ! note that this re-sets m%udotdot_TP and m%UFL, but they are the same values as earlier in this routine so it doesn't change results in SDOut_MapOutputs()
            CALL SD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat2, ErrMsg2 ); if(Failed()) return
            !Assign the acceleration to the x variable since it will be used for output file purposes for SSqmdd01-99, and dxdt will disappear
            m%qmdotdot=dxdt%qmdot
            ! Destroy dxdt because it is not necessary for the rest of the subroutine
            CALL SD_DestroyContState( dxdt, ErrStat2, ErrMsg2); if(Failed()) return
         END IF
          
         ! Write the previous output data into the output file           
         IF ( ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) .AND. ( t > m%LastOutTime ) ) THEN
            IF ((m%Decimat .EQ. p%OutDec) .OR. (m%Decimat .EQ. 0))  THEN
               m%Decimat=1  !reset counter
               CALL SDOut_WriteOutputs( p%UnJckF, m%LastOutTime, m%SDWrOutput, p, ErrStat2, ErrMsg2 ); if(Failed()) return
            ELSE      
               m%Decimat=m%Decimat+1
            ENDIF
         END IF        
         
         ! Map calculated results into the AllOuts Array + perform averaging and all necessary extra calculations
         CALL SDOut_MapOutputs(t, u,p,x,y, m, AllOuts, ErrStat2, ErrMsg2); if(Failed()) return
            
         ! Put the output data in the WriteOutput array
         DO I = 1,p%NumOuts+p%OutAllInt*p%OutAllDims
            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
            IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN
               m%SDWrOutput(I) = y%WriteOutput(I)            
            END IF                        
         END DO
         m%LastOutTime   = t
      ENDIF           
  
CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed
   
   SUBROUTINE CleanUp
       CALL SD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
   END SUBROUTINE CleanUp

END SUBROUTINE SD_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
!! note that this also sets m%UFL and m%udotdot_TP
SUBROUTINE SD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
      REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(SD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
      TYPE(SD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType), INTENT(IN)     :: x           !< Continuous states at t
      TYPE(SD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(SD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(SD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      TYPE(SD_ContinuousStateType), INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at t
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      REAL(ReKi)                                   :: rotations(3)
      INTEGER(IntKi)       :: ErrStat2
      CHARACTER(ErrMsgLen) :: ErrMsg2
      REAL(ReKi), ALLOCATABLE                      :: seismic_force_temp(:)
      CHARACTER(*), PARAMETER                      :: RoutineName = 'SD_CalcContStateDeriv'

      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      IF (p%CBMod) THEN  !Craig-Bampton reduction
          
         ! INTENT(OUT) automatically deallocates the arrays on entry, we have to allocate them here
         CALL AllocAry(dxdt%qm,    p%qmL, 'dxdt%qm',    ErrStat2, ErrMsg2 ); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcContStateDeriv' )
         CALL AllocAry(dxdt%qmdot, p%qmL, 'dxdt%qmdot', ErrStat2, ErrMsg2 ); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcContStateDeriv' )

         IF ( ErrStat >= AbortErrLev ) RETURN
         
         IF ( p%qmL == 0 ) RETURN

         ! form u(2)
         m%udot_TP    = (/u%TPMesh%TranslationVel( :,1), u%TPMesh%RotationVel(:,1)/)
      
         ! form u(3) in Eq. 10:
         m%udotdot_TP = (/u%TPMesh%TranslationAcc(:,1), u%TPMesh%RotationAcc(:,1)/)
      
         ! form u(4) in Eq. 10:
         CALL ConstructUFL( u, p, m%UFL )
      
         !Equation 12: X=A*x + B*u + Fx (Eq 12)
         dxdt%qm= x%qmdot

         ! NOTE: matmul( TRANSPOSE(p%PhiM), m%UFL ) = matmul( m%UFL, p%PhiM ) because UFL is 1-D
                != a(2,1) * x(1)   +   a(2,2) * x(2)         +    b(2,2) * u(2)            +  b(2,3) * u(3)                       + b(2,4) * u(4)             + fx(2) 

         dxdt%qmdot = p%NOmegaM2*x%qm                                   &
                    + p%N2OmegaMJDamp*x%qmdot - matmul(p%CMM , x%qmdot) &
                    - matmul(p%CMB , m%udot_TP)                         &
                    - matmul(p%MMB , m%udotdot_TP)                      &
                    + matmul(m%UFL , p%PhiM ) + p%FX 

         IF (p%SeismicInp) THEN
           dxdt%qmdot = dxdt%qmdot + m%Ug * p%FMIMKP + m%Udotg * p%FMIMCP - m%Uddotg * p%FMIMMP + m%PHIg * p%FMIMKP_PHI + m%PHIdotg * p%FMIMCP_PHI - m%PHIddotg * p%FMIMMP_PHI + m%Vg * p%FMIMKP_V + m%Vdotg * p%FMIMCP_V - m%Vddotg * p%FMIMMP_V
         ENDIF

      ELSE  !full fem without modal reduction

         CALL AllocAry(seismic_force_temp, p%DOFL, 'seismic_force_temp', ErrStat2, ErrMsg2 ); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcContStateDeriv' )

         ! Compute the small rotation angles given the input direction cosine matrix
         rotations  = GetSmllRotAngs(u%TPMesh%Orientation(:,:,1), ErrStat2, Errmsg2); if(Failed()) return

         ! form u(1)
         m%u_TP = (/REAL(u%TPMesh%TranslationDisp(:,1),ReKi), rotations/)

         ! form u(2)
         m%udot_TP    = (/u%TPMesh%TranslationVel( :,1), u%TPMesh%RotationVel(:,1)/)
      
         ! form u(3)
         m%udotdot_TP = (/u%TPMesh%TranslationAcc(:,1), u%TPMesh%RotationAcc(:,1)/)
      
         ! form u(4)
         CALL ConstructUFL( u, p, m%UFL )

         !Equation 1 in X=A*x + B*u + Fx
         dxdt%qm= x%qmdot

         !Equation 2 in X=A*x + B*u + Fx
         dxdt%qmdot = matmul(p%A_21 , x%qm      )    &
                    + matmul(p%A_22 , x%qmdot   )    &
                    + matmul(p%B_21 , m%u_TP    )    &
                    + matmul(p%B_22 , m%udot_TP )    &
                    + matmul(p%B_23 , m%udotdot_TP ) &
                    + matmul(p%B_24 , m%UFL        ) &
                    + p%FX                           

          IF (p%SeismicInp) THEN

             seismic_force_temp = m%Ug * matmul(p%B_24 , p%FSISKL)
             dxdt%qmdot = dxdt%qmdot + seismic_force_temp

             seismic_force_temp = m%Udotg * matmul(p%B_24 , p%FSISCL)
             dxdt%qmdot = dxdt%qmdot + seismic_force_temp

             seismic_force_temp = m%Uddotg * matmul(p%B_24 , p%FSISML)
             dxdt%qmdot = dxdt%qmdot - seismic_force_temp
 
         ENDIF      

         DEALLOCATE(seismic_force_temp)  

      ENDIF !Craig-Bampton vs full fem

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed

END SUBROUTINE SD_CalcContStateDeriv

!-----------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_Input(SDInputFile, Init, p, ErrStat,ErrMsg)
   CHARACTER(*),            INTENT(IN)     :: SDInputFile
   TYPE(SD_InitType) ,      INTENT(INOUT)  :: Init
   TYPE(SD_ParameterType) , INTENT(INOUT)  :: p
   INTEGER(IntKi),          INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),            INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None
! local variable for input and output
CHARACTER(1024)              :: PriPath                                         ! The path to the primary input file
CHARACTER(1024)              :: Line                                            ! String to temporarially hold value of read line
INTEGER                      :: Sttus

LOGICAL                      :: Echo  
INTEGER(IntKi)               :: UnIn, UnInUg
INTEGER(IntKi)               :: IOS
INTEGER(IntKi)               :: UnEc   !Echo file ID

REAL(ReKi),PARAMETER        :: WrongNo=-9999.   ! Placeholder value for bad(old) values in JDampings

INTEGER(IntKi)               :: I, J, flg, K
REAL(ReKi)                   :: Dummy_ReAry(SDMaxInpCols) 
INTEGER(IntKi)               :: Dummy_IntAry(SDMaxInpCols)
INTEGER(IntKi)               :: ErrStat2
CHARACTER(ErrMsgLen)         :: ErrMsg2
CHARACTER(1024)              :: UgFile        !  File that contains the seismic input data
! Initialize ErrStat
ErrStat = ErrID_None
ErrMsg  = ""

UnEc = -1 
Echo = .FALSE.

CALL GetNewUnit( UnIn )   
  
CALL OpenFInpfile(UnIn, TRIM(SDInputFile), ErrStat2, ErrMsg2)

IF ( ErrStat2 /= ErrID_None ) THEN
   Call Fatal('Could not open SubDyn input file')
   return
END IF

CALL GetPath( SDInputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.


!-------------------------- HEADER ---------------------------------------------
CALL ReadCom( UnIn, SDInputFile, 'SubDyn input file header line 1', ErrStat2, ErrMsg2 ); if(Failed()) return
CALL ReadCom( UnIn, SDInputFile, 'SubDyn input file header line 2', ErrStat2, ErrMsg2 ); if(Failed()) return

!-------------------------- SIMULATION CONTROL PARAMETERS ----------------------
CALL ReadCom( UnIn, SDInputFile, ' SIMULATION CONTROL PARAMETERS ', ErrStat2, ErrMsg2 ); if(Failed()) return
CALL ReadVar(UnIn, SDInputFile, Echo, 'Echo', 'Echo Input File Logic Variable',ErrStat2, ErrMsg2); if(Failed()) return

IF ( Echo )  THEN 
   CALL OpenEcho ( UnEc, TRIM(Init%RootName)//'.ech' ,ErrStat2, ErrMsg2)
   IF ( ErrStat2 /= 0 ) THEN
      CALL Fatal("Could not open SubDyn echo file")
      return
   END IF
   REWIND(UnIn)
   !bjj: note we don't need to do error checking here; it was already checked (this is just a repeat of above)
   CALL ReadCom( UnIn, SDInputFile, 'SubDyn input file header line 1', ErrStat2, ErrMsg2 )
   CALL ReadCom( UnIn, SDInputFile, 'SubDyn input file header line 2', ErrStat2, ErrMsg2 )
   CALL ReadCom( UnIn, SDInputFile, 'SIMULATION CONTROL PARAMETERS'  , ErrStat2, ErrMsg2, UnEc )
   CALL ReadVar( UnIn, SDInputFile, Echo, 'Echo', 'Echo Input File Logic Variable',ErrStat2, ErrMsg2, UnEc )
ENDIF 

! Read time step   ("default" means use the glue-code default)
CALL ReadVar( UnIn, SDInputFile, Line, 'SDdeltaT', 'Subdyn Time Step',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

CALL Conv2UC( Line )    ! Convert Line to upper case.
IF ( TRIM(Line) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value timestep provided by the glue code.
    p%SDdeltaT=Init%DT
ELSE                                   ! The input must have been specified numerically.
   READ (Line,*,IOSTAT=IOS)  p%SDdeltaT
   CALL CheckIOS ( IOS, SDInputFile, 'SDdeltaT', NumType, ErrStat2,ErrMsg2 ); if(Failed()) return

   IF ( ( p%SDdeltaT <=  0 ) )  THEN 
      call Fatal('SDdeltaT must be greater than or equal to 0.')
      return         
   END IF  
END IF
      
CALL ReadVar ( UnIn, SDInputFile, p%IntMethod, 'IntMethod', 'Integration Method',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadLVar(UnIn, SDInputFile, p%SttcSolve, 'SttcSolve', 'Solve dynamics about static equilibrium point', ErrStat2, ErrMsg2, UnEc); if(Failed()) return
CALL ReadLVar(UnIn, SDInputFile, p%SeismicInp, 'SeismicInp', 'Existence of seismic input Ug', ErrStat2, ErrMsg2, UnEc); if(Failed()) return
!-------------------- FEA and CRAIG-BAMPTON PARAMETERS---------------------------
CALL ReadCom  ( UnIn, SDInputFile, ' FEA and CRAIG-BAMPTON PARAMETERS ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%FEMMod, 'FEMMod', 'FEM analysis mode'             ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return ! 0= Euler-Bernoulli(E-B); 1=Tapered E-B; 2= Timoshenko; 3= tapered Timoshenko
CALL ReadIVar ( UnIn, SDInputFile, Init%NDiv  , 'NDiv'  , 'Number of divisions per member',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadLVar ( UnIn, SDInputFile, Init%CBMod , 'CBMod' , 'C-B mod flag'                  ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
p%CBmod = Init%CBMod !We copy this variable because Init is not available in SD_CalcOutput y SD_CalcContStateDeriv, and we now need this information
                     !better than adding new arguments to those routines, that are called from different points.

IF (Check( (p%IntMethod < 1) .OR.(p%IntMethod > 4)     , 'IntMethod must be 1 through 4.')) return
IF (Check( (Init%FEMMod < 0 ) .OR. ( Init%FEMMod > 4 ) , 'FEMMod must be 0, 1, 2, or 3.')) return
IF (Check( Init%NDiv < 1                               , 'NDiv must be a positive integer')) return

IF (Init%CBMod) THEN
   ! Nmodes - Number of interal modes to retain.
   CALL ReadIVar ( UnIn, SDInputFile, p%Nmodes, 'Nmodes', 'Number of internal modes',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

   IF (Check( p%Nmodes < 0 , 'Nmodes must be a non-negative integer.')) return
   
   if ( p%Nmodes > 0 ) THEN
      ! Damping ratios for retained modes
      CALL AllocAry(Init%JDampings, p%Nmodes, 'JDamping', ErrStat2, ErrMsg2) ; if(Failed()) return
      Init%JDampings=WrongNo !Initialize
   
      CALL ReadAry( UnIn, SDInputFile, Init%JDampings, p%Nmodes, 'JDamping', 'Damping ratio of the internal modes', ErrStat2, ErrMsg2, UnEc );
      ! note that we don't check the ErrStat2 here; if the user entered fewer than Nmodes values, we will use the
      ! last entry to fill in remaining values.
      !Check 1st value, we need at least one good value from user or throw error
      IF ((Init%JDampings(1) < 0 ) .OR. (Init%JDampings(1) >= 100.0)) THEN
            CALL Fatal('Damping ratio should be larger than 0 and less than 100')
            return
      ELSE
         DO I = 2, p%Nmodes
            IF ( Init%JDampings(I) .EQ. WrongNo ) THEN
               Init%Jdampings(I:p%Nmodes)=Init%JDampings(I-1)
               IF (i /= 2) THEN ! display an informational message if we're repeating the last value (unless we only entered one value)
                  ErrStat = ErrID_Info
                  ErrMsg  = 'Using damping ratio '//trim(num2lstr(Init%JDampings(I-1)))//' for modes '//trim(num2lstr(I))//' - '//trim(num2lstr(p%Nmodes))//'.'
               END IF
               EXIT
            ELSEIF ( ( Init%JDampings(I) < 0 ) .OR.( Init%JDampings(I) >= 100.0 ) ) THEN    
               CALL Fatal('Damping ratio should be larger than 0 and less than 100')
               return
            ENDIF      
        ENDDO
      ENDIF   
      IF (ErrStat2 /= ErrID_None .AND. Echo) THEN ! ReadAry had an error because it couldn't read the entire array so it didn't write this to the echo file; we assume the last-read values are used for remaining JDampings
         WRITE( UnEc, Ec_ReAryFrmt ) 'JDamping', 'Damping ratio of the internal modes', Init%Jdampings(1:MIN(p%Nmodes,NWTC_MaxAryLen))              
      END IF
   ELSE
      CALL ReadCom( UnIn, SDInputFile, 'JDamping', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   END IF

ELSE   !CBMOD=FALSE  : all modes are retained, not sure how many they are yet
   !note at this stage I do not know DOFL yet; Nmodes will be updated later for the FULL FEM CASE. 
   p%Nmodes = -1
   !Ignore next line
   CALL ReadCom( UnIn, SDInputFile, 'Nmodes', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   !Read 1 damping value for all modes
   CALL AllocAry(Init%JDampings, 1, 'JDamping', ErrStat2, ErrMsg2) ; if(Failed()) return
   CALL ReadVar ( UnIn, SDInputFile, Init%JDampings(1), 'JDampings', 'Damping ratio',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   IF ( ( Init%JDampings(1) < 0 ) .OR.( Init%JDampings(1) >= 100.0 ) ) THEN 
         CALL Fatal('Damping ratio should be larger than 0 and less than 100.')
         RETURN
   ENDIF
ENDIF

IF ((p%Nmodes > 0) .OR. (.NOT.(Init%CBMod))) THEN !This if should not be at all, dampings should be divided by 100 regardless, also if CBmod=false p%Nmodes is undefined, but if Nmodes=0 then JDampings does not exist
   Init%JDampings = Init%JDampings/100.0_ReKi   !now the 20 is .20 as it should in all cases for 1 or Nmodes JDampings
END IF

!--------------------- STRUCTURE JOINTS: joints connect structure members -------------------------------
CALL ReadCom  ( UnIn, SDInputFile,               'STRUCTURE JOINTS'           ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NJoints, 'NJoints', 'Number of joints',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,               'Joint Coordinates Headers'  ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,               'Joint Coordinates Units'    ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%Joints, Init%NJoints, JointsCol, 'Joints', ErrStat2, ErrMsg2 ); if(Failed()) return
DO I = 1, Init%NJoints
   CALL ReadAry( UnIn, SDInputFile, Dummy_ReAry, JointsCol, 'Joints', 'Joint number and coordinates', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   Init%Joints(I,:) = Dummy_ReAry(1:JointsCol)
ENDDO
IF (Check(  Init%NJoints < 2, 'NJoints must be greater than 1')) return

!---------- GO AHEAD  and ROTATE STRUCTURE IF DESIRED TO SIMULATE WINDS FROM OTHER DIRECTIONS -------------
CALL SubRotate(Init%Joints,Init%NJoints,Init%SubRotateZ)

!------------------- BASE REACTION JOINTS: T/F for Locked/Free DOF @ each Reaction Node ---------------------
! The joints should be all clamped for now 
CALL ReadCom  ( UnIn, SDInputFile,           'BASE REACTION JOINTS'                           ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, p%NReact, 'NReact', 'Number of joints with reaction forces',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,           'Base reaction joints headers '                  ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,           'Base reaction joints units   '                  ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(p%Reacts, p%NReact, ReactCol, 'Reacts', ErrStat2, ErrMsg2 ); if(Failed()) return
DO I = 1, p%NReact
   CALL ReadAry( UnIn, SDInputFile, Dummy_IntAry, ReactCol, 'Reacts', 'Joint number and dof', ErrStat2 ,ErrMsg2, UnEc); if(Failed()) return
   p%Reacts(I,:) = Dummy_IntAry(1:ReactCol)
ENDDO
IF (Check ( ( p%NReact < 1 ) .OR. (p%NReact > Init%NJoints) , 'NReact must be greater than 0 and less than number of joints')) return

!------- INTERFACE JOINTS: T/F for Locked (to the TP)/Free DOF @each Interface Joint (only Locked-to-TP implemented thus far (=rigid TP)) ---------
! Joints with reaction forces, joint number and locked/free dof
CALL ReadCom  ( UnIn, SDInputFile,               'INTERFACE JOINTS'                     ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NInterf, 'NInterf', 'Number of joints fixed to TP',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,               'Interface joints headers',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,               'Interface joints units  ',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%Interf, Init%NInterf, InterfCol, 'Interf', ErrStat2, ErrMsg2); if(Failed()) return
DO I = 1, Init%NInterf
   CALL ReadIAry( UnIn, SDInputFile, Dummy_IntAry, InterfCol, 'Interf', 'Interface joint number and dof', ErrStat2,ErrMsg2, UnEc); if(Failed()) return
   Init%Interf(I,:) = Dummy_IntAry(1:InterfCol)
ENDDO
IF (Check( ( Init%NInterf < 0 ) .OR. (Init%NInterf > Init%NJoints), 'NInterf must be non-negative and less than number of joints.')) RETURN

!----------------------------------- MEMBERS --------------------------------------
! One day we will need to take care of COSMIDs for non-circular members
CALL ReadCom  ( UnIn, SDInputFile,             'Members '                     ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, p%NMembers, 'NMembers', 'Number of members',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,             'Members Headers'              ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,             'Members Units  '              ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%Members, p%NMembers, MembersCol, 'Members', ErrStat2, ErrMsg2)
p%NSLPMEl = 0
DO I = 1, p%NMembers
   CALL ReadAry( UnIn, SDInputFile, Dummy_IntAry, MembersCol, 'Members', 'Member number and connectivity ', ErrStat2,ErrMsg2, UnEc); if(Failed()) return
   Init%Members(I,:) = Dummy_IntAry(1:MembersCol)
   p%NSLPMEl = p%NSLPMEl + Init%Members(I,6) ! Count number of SLPM Elements
ENDDO   
IF (Check( p%NMembers < 1 , 'NMembers must be > 0')) return

!------------------ MEMBER X-SECTION PROPERTY data 1/2 [isotropic material for now: use this table if circular-tubular elements ------------------------
CALL ReadCom  ( UnIn, SDInputFile,                 ' Member X-Section Property Data 1/2 ',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NPropSets, 'NPropSets', 'Number of property sets',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                 'Property Data 1/2 Header'            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                 'Property Data 1/2 Units '            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%PropSets, Init%NPropSets, PropSetsCol, 'ProSets', ErrStat2, ErrMsg2) ; if(Failed()) return
DO I = 1, Init%NPropSets
   CALL ReadAry( UnIn, SDInputFile, Dummy_ReAry, PropSetsCol, 'PropSets', 'PropSets number and values ', ErrStat2 , ErrMsg2, UnEc); if(Failed()) return
   Init%PropSets(I,:) = Dummy_ReAry(1:PropSetsCol)
ENDDO   
IF (Check( Init%NPropSets < 1 , 'NPropSets must be >0')) return

!------------------ MEMBER X-SECTION PROPERTY data 2/2 [isotropic material for now: use this table if any section other than circular, however provide COSM(i,j) below) ------------------------
CALL ReadCom  ( UnIn, SDInputFile,                  'Member X-Section Property Data 2/2 '               ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NXPropSets, 'NXPropSets', 'Number of non-circular property sets',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                  'Property Data 2/2 Header'                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                  'Property Data 2/2 Unit  '                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%XPropSets, Init%NXPropSets, XPropSetsCol, 'XPropSets', ErrStat2, ErrMsg2); if(Failed()) return
DO I = 1, Init%NXPropSets
   CALL ReadAry( UnIn, SDInputFile, Init%XPropSets(I,:), XPropSetsCol, 'XPropSets', 'XPropSets ID and values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
ENDDO   
IF (Check( Init%NXPropSets < 0, 'NXPropSets must be >=0')) return

!------------------ SIMPLIFIED LPM-TYPE MEMBERS PROPERTIES  ------------------------
CALL ReadCom  ( UnIn, SDInputFile,                  'Simplified LPM-Tye Members Property Data'          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NSLPMPropSets, 'SLPMPropSets','Number of Simplified LPM property sets',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                  'SLPM Property Data Header'                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                  'SLPM Property Data Unit  '                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%SLPMPropSets, Init%NSLPMPropSets, SLPMPropSetsCol, 'SLPMPropSets', ErrStat2, ErrMsg2); if(Failed()) return
DO I = 1, Init%NSLPMPropSets
   CALL ReadAry( UnIn, SDInputFile, Init%SLPMPropSets(I,:), SLPMPropSetsCol, 'SLPMPropSets', 'SLPMPropSets ID and values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
ENDDO   
IF (Check( Init%NSLPMPropSets < 0, 'NSLPMPropSets must be >=0')) return

!---------------------- MEMBER COSINE MATRICES COSM(i,j) ------------------------
CALL ReadCom  ( UnIn, SDInputFile,              'Member direction cosine matrices '                   ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NCOSMs, 'NCOSMs', 'Number of unique direction cosine matrices',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,              'Cosine Matrices Headers'                             ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,              'Cosine Matrices Units  '                             ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%COSMs, Init%NCOSMs, COSMsCol, 'COSMs', ErrStat2, ErrMsg2); if(Failed()) return
DO I = 1, Init%NCOSMs
   CALL ReadAry( UnIn, SDInputFile, Init%COSMs(I,:), COSMsCol, 'CosM', 'Cosine Matrix IDs  and Values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
ENDDO   
IF (Check( Init%NCOSMs < 0     ,'NCOSMs must be >=0')) return

!------------------------ JOINT ADDITIONAL CONCENTRATED MASSES--------------------------
CALL ReadCom  ( UnIn, SDInputFile,              'Additional concentrated masses at joints '               ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NCMass, 'NCMass', 'Number of joints that have concentrated masses',ErrStat2, ErrMsg2, UnEc); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,              'Concentrated Mass Headers'                               ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,              'Concentrated Mass Units'                                 ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%CMass, Init%NCMass, CMassCol, 'CMass', ErrStat2, ErrMsg2); if(Failed()) return
Init%CMass = 0.0
DO I = 1, Init%NCMass
   CALL ReadAry( UnIn, SDInputFile, Init%CMass(I,:), CMassCol, 'CMass', 'Joint number and mass values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
ENDDO   
IF (Check( Init%NCMass < 0     , 'NCMass must be >=0')) return

!------------------------ SEISMIC INPUT --------------------------
CALL ReadCom  ( UnIn, SDInputFile,              'Seismic Input'               ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
IF (p%SeismicInp) THEN
      ! UgFile - Name of the file containing earthquake signal
   CALL ReadVar ( UnIn, SDInputFile, UgFile, 'UgFile', 'Name of the file containing earthquake signal', ErrStat2, ErrMsg2, UnEc )
   IF ( PathIsRelative( UgFile ) ) UgFile = TRIM(PriPath)//TRIM(UgFile)
   CALL GetNewUnit( UnInUg )   
   CALL OpenFInpfile(UnInUg, TRIM(UgFile), ErrStat2, ErrMsg2)
   IF ( ErrStat2 /= ErrID_None ) THEN
      Call Fatal('Could not open SubDyn seismic signal input file')
      return
   END IF
       !-------------------------- HEADER ---------------------------------------------
   CALL ReadCom( UnInUg, UgFile, 'SubDyn seismic signal input file header line 1', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadCom( UnInUg, UgFile, 'SubDyn seismic signal input file header line 2', ErrStat2, ErrMsg2 ); if(Failed()) return
       !-------------------------- SEISMIC SIGNAL ----------------------
   CALL ReadCom( UnInUg, UgFile, ' SEISMIC SIGNAL ', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVar( UnInUg, UgFile, p%SDDeltaTUg, 'SDdeltaTUg', 'Seismic signal Time Step',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   CALL ReadIVar (UnInUg, UgFile, p%NtUg, 'NtUg', 'Number of data points',ErrStat2, ErrMsg2, UnEc); if(Failed()) return
   CALL ReadVar( UnInUg, UgFile, p%UgDir, 'UgDir', 'Shaking direction',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   CALL AllocAry(p%UgData, p%NtUg, 10, 'UgData', ErrStat2, ErrMsg2); if(Failed()) return
   CALL ReadCom( UnInUg, UgFile, 'SubDyn seismic signal input file header line 3', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadCom( UnInUg, UgFile, 'units', ErrStat2, ErrMsg2 ); if(Failed()) return
   DO I = 1, p%NtUg
      CALL ReadAry( UnInUg, UgFile, p%UgData(I,:), 10, 'UgData', 'time, input ground displacement, velocity and acceleration, input ground twist, angular velocity and acceleration, vertical displacement, velocity and acceleration', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   ENDDO
   CLOSE( UnInUg )
ENDIF
!---------------------------- OUTPUT: SUMMARY & OUTFILE ------------------------------
CALL ReadCom (UnIn, SDInputFile,               'OUTPUT'                                            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadLVar(UnIn, SDInputFile, Init%SSSum  , 'SSSum'  , 'Summary File Logic Variable'            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadLVar(UnIn, SDInputFile, Init%OutCOSM, 'OutCOSM', 'Cosine Matrix Logic Variable'           ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return !bjj: TODO: OutCOSM isn't used anywhere else.
CALL ReadLVar(UnIn, SDInputFile, p%OutAll    , 'OutAll' , 'Output all Member Forces Logic Variable',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
!Store an integer version of it
p%OutAllInt= 1
IF ( .NOT. p%OutAll ) p%OutAllInt= 0
CALL ReadIVar(UnIn, SDInputFile, p%OutSwtch, 'OutSwtch', 'Output to which file variable',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
IF (Check( ( p%OutSwtch < 1 ) .OR. ( p%OutSwtch > 3) ,'OutSwtch must be >0 and <4')) return

Swtch: SELECT CASE (p%OutSwtch)
 CASE (1, 3) Swtch
    !p%OutJckF = TRIM(Init%RootName)//'.out'
 CASE (2)  Swtch
    !pass to glue code
 CASE DEFAULT Swtch
    CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": OutSwtch must be >0 and <4')
    return
 END SELECT Swtch
     
! TabDelim - Output format for tabular data.
CALL ReadLVar ( UnIn,  SDInputFile, Init%TabDelim, 'TabDelim', 'Use Tab Delimitation for numerical outputs',ErrStat2, ErrMsg2, UnEc); if(Failed()) return
IF ( Init%TabDelim ) THEN
         p%Delim = TAB
ELSE
         p%Delim = ' '
END IF

CALL ReadIVar( UnIn, SDInputFile, p%OutDec  , 'OutDec'  , 'Output Decimation'                , ErrStat2 , ErrMsg2 , UnEc ); if(Failed()) return
CALL ReadVar ( UnIn, SDInputFile, p%OutFmt  , 'OutFmt'  , 'Format for numerical outputs'     , ErrStat2 , ErrMsg2 , UnEc ); if(Failed()) return
CALL ReadVar ( UnIn, SDInputFile, p%OutSFmt , 'OutSFmt' , 'Format for output column headers' , ErrStat2 , ErrMsg2 , UnEc ); if(Failed()) return
CALL ReadCom ( UnIn, SDInputFile,             ' Member Output List SECTION ',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar( UnIn, SDInputFile, p%NMOutputs, 'NMOutputs', 'Number of Members whose output must go into OutJckF and/or FAST .out',ErrStat2, ErrMsg2, UnEc )
if (Failed()) return
IF (Check ( (p%NMOutputs < 0) .OR. (p%NMOutputs > p%NMembers) .OR. (p%NMOutputs > 9), 'NMOutputs must be >=0 and <= minimim(NMembers,9)')) return

CALL ReadCom( UnIn, SDInputFile, ' Output Member Headers',ErrStat2, ErrMsg2, UnEc) ; if(Failed()) return
CALL ReadCom( UnIn, SDInputFile, ' Output Member Units'  ,ErrStat2, ErrMsg2, UnEc) ; if(Failed()) return

IF ( p%NMOutputs > 0 ) THEN
   ! Allocate memory for filled group arrays
   ALLOCATE ( p%MOutLst(p%NMOutputs), STAT = ErrStat2 )     !this list contains different arrays for each of its elements
   IF ( ErrStat2 /= ErrID_None ) THEN
      CALL  Fatal(' Error in file "'//TRIM(SDInputFile)//': Error allocating MOutLst arrays')
      RETURN
   END IF

   DO I = 1,p%NMOutputs
      READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line 
      IF (ErrStat2 == 0) THEN
         READ(Line,*,IOSTAT=ErrStat2) p%MOutLst(I)%MemberID, p%MOutLst(I)%NOutCnt
         IF ( ErrStat2 /= 0 .OR. p%MOutLst(I)%NOutCnt < 1 .OR. p%MOutLst(I)%NOutCnt > 9 .OR. p%MOutLst(I)%NOutCnt > Init%Ndiv+1) THEN
            CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": NOutCnt must be >= 1 and <= minimim(Ndiv+1,9)')
            RETURN
         END IF            
         CALL AllocAry( p%MOutLst(I)%NodeCnt, p%MOutLst(I)%NOutCnt, 'NodeCnt', ErrStat2, ErrMsg2); if(Failed()) return

         READ(Line,*,IOSTAT=ErrStat2) p%MOutLst(I)%MemberID,  p%MOutLst(I)%NOutCnt,  p%MOutLst(I)%NodeCnt
         IF ( Check( ErrStat2 /= 0 , 'Failed to read member output list properties.')) return

         ! Check if MemberID is in the member list and the NodeCnt is a valid number
         flg = 0
         DO J = 1, p%NMembers
            IF(p%MOutLst(I)%MemberID .EQ. Init%Members(j, 1)) THEN
               flg = flg + 1 ! flg could be greater than 1, when there are more than 9 internal nodes of a member.
               IF( (p%MOutLst(I)%NOutCnt < 10) .and. ((p%MOutLst(I)%NOutCnt > 0)) ) THEN
                  DO K = 1,p%MOutLst(I)%NOutCnt
                     ! node number should be less than NDiv + 1
                     IF( (p%MOutLst(I)%NodeCnt(k) > (Init%NDiv+1)) .or. (p%MOutLst(I)%NodeCnt(k) < 1) ) THEN
                        CALL Fatal(' NodeCnt should be less than NDIV+1 and greater than 0. ')
                        RETURN
                     ENDIF
                  ENDDO
               ELSE
                  CALL Fatal(' NOutCnt should be less than 10 and greater than 0. ')
                  RETURN
               ENDIF
            ENDIF
         ENDDO
         IF (Check (flg .EQ. 0 , ' MemberID is not in the Members list. ')) return

         IF ( Echo ) THEN
            WRITE( UnEc, '(A)' ) TRIM(Line)
         END IF
      END IF
   END DO
END IF 

! OutList - list of requested parameters to output to a file
CALL ReadCom( UnIn, SDInputFile, 'SSOutList',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

ALLOCATE(Init%SSOutList(MaxOutChs), STAT=ErrStat2)
If (Check( ErrStat2 /= ErrID_None ,'Error allocating SSOutList arrays')) return

CALL ReadOutputList ( UnIn, SDInputFile, Init%SSOutList, p%NumOuts, &
                                              'SSOutList', 'List of outputs requested', ErrStat2, ErrMsg2, UnEc )
if(Failed()) return
CALL CleanUp()

CONTAINS

   LOGICAL FUNCTION Check(Condition, ErrMsg_in)
        logical, intent(in) :: Condition
        character(len=*), intent(in) :: ErrMsg_in
        Check=Condition
        if (Check) call Fatal(' Error in file '//TRIM(SDInputFile)//': '//trim(ErrMsg_in))
   END FUNCTION Check

   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Input') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed

   SUBROUTINE Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'SD_Input');
      CALL CleanUp()
   END SUBROUTINE Fatal

   SUBROUTINE CleanUp()
      CLOSE( UnIn )
      IF (Echo) CLOSE( UnEc )
   END SUBROUTINE

END SUBROUTINE SD_Input


!----------------------------------------------------------------------------------------------------------------------------------
!> Rotate the joint coordinates with respect to global z
SUBROUTINE SubRotate(Joints,NJoints,SubRotZ)
   REAL(ReKi),                       INTENT(IN)       :: SubRotZ    ! Rotational angle in degrees
   INTEGER(IntKi),                   INTENT(IN)       :: NJOINTS    ! Row size of Joints 
   REAL(ReKi), DIMENSION(NJOINTS,3), INTENT(INOUT)    :: JOINTS     ! Rotational angle in degrees (Njoints,4)
   !locals
   REAL(ReKi)                 :: rot  !angle in rad
   REAL(ReKi), DIMENSION(2,2) :: ROTM !rotational matrix (cos matrix with -theta)
   
   rot=pi*SubRotz/180.
   ROTM=transpose(reshape([ COS(rot),    -SIN(rot) , &
                            SIN(rot) ,    COS(rot)], [2,2] ))
   Joints(:,2:3)= transpose(matmul(ROTM,transpose(Joints(:,2:3))))

END SUBROUTINE  SubRotate           

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE SD_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      TYPE(SD_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(SD_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(SD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(SD_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(SD_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(SD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states            
      TYPE(SD_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(SD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ! Initialize ErrStat
      ErrStat = ErrID_None         
      ErrMsg  = ""               

      ! Determine if we need to close the output file
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN   
         IF ((m%Decimat .EQ. p%OutDec) .OR. (m%Decimat .EQ. 0))  THEN
               ! Write out the last stored set of outputs before closing
            CALL SDOut_WriteOutputs( p%UnJckF, m%LastOutTime, m%SDWrOutput, p, ErrStat, ErrMsg )   
         ENDIF
         CALL SDOut_CloseOutput( p, ErrStat, ErrMsg )         
      END IF 
      
      ! Destroy data
      CALL SD_DestroyInput( u, ErrStat, ErrMsg )
      CALL SD_DestroyParam( p, ErrStat, ErrMsg )
      CALL SD_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL SD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL SD_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL SD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
      CALL SD_DestroyMisc( m,  ErrStat, ErrMsg )
      CALL SD_DestroyOutput( y, ErrStat, ErrMsg )

END SUBROUTINE SD_End

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
!! equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!  See, e.g.,
!!    - http://en.wikipedia.org/wiki/Linear_multistep_method
!!    - K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE SD_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other states at t on input at t + dt on output
      TYPE(SD_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ! local variables
      TYPE(SD_ContinuousStateType) :: xdot       ! Continuous state derivs at t
      TYPE(SD_InputType)           :: u_interp

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! need xdot at t
      CALL SD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg  )  ! we need to allocate input arrays/meshes before calling ExtrapInterp...
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
      CALL SD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat, ErrMsg ) ! initializes xdot
      CALL SD_DestroyInput( u_interp, ErrStat, ErrMsg)   ! we don't need this local copy anymore

      if (n <= 2) then
         OtherState%n = n
         !OtherState%xdot ( 3 - n ) = xdot
         CALL SD_CopyContState( xdot, OtherState%xdot ( 3 - n ), MESH_UPDATECOPY, ErrStat, ErrMsg )
         CALL SD_RK4(t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      else
         if (OtherState%n < n) then
            OtherState%n = n
            CALL SD_CopyContState( OtherState%xdot ( 3 ), OtherState%xdot ( 4 ), MESH_UPDATECOPY, ErrStat, ErrMsg )
            CALL SD_CopyContState( OtherState%xdot ( 2 ), OtherState%xdot ( 3 ), MESH_UPDATECOPY, ErrStat, ErrMsg )
            CALL SD_CopyContState( OtherState%xdot ( 1 ), OtherState%xdot ( 2 ), MESH_UPDATECOPY, ErrStat, ErrMsg )
            !OtherState%xdot(4)    = OtherState%xdot(3)
            !OtherState%xdot(3)    = OtherState%xdot(2)
            !OtherState%xdot(2)    = OtherState%xdot(1)
         elseif (OtherState%n > n) then
            ErrStat = ErrID_Fatal
            ErrMsg = ' Backing up in time is not supported with a multistep method '
            RETURN
         endif
         CALL SD_CopyContState( xdot, OtherState%xdot ( 1 ), MESH_UPDATECOPY, ErrStat, ErrMsg )
         !OtherState%xdot ( 1 )     = xdot  ! make sure this is most up to date
         x%qm    = x%qm    + (p%SDDeltaT / 24.) * ( 55.*OtherState%xdot(1)%qm - 59.*OtherState%xdot(2)%qm    + 37.*OtherState%xdot(3)%qm  &
                                       - 9. * OtherState%xdot(4)%qm )
         x%qmdot = x%qmdot + (p%SDDeltaT / 24.) * ( 55.*OtherState%xdot(1)%qmdot - 59.*OtherState%xdot(2)%qmdot  &
                                          + 37.*OtherState%xdot(3)%qmdot  - 9.*OtherState%xdot(4)%qmdot )
      endif
      CALL SD_DestroyContState(xdot, ErrStat, ErrMsg)
      CALL SD_DestroyInput(u_interp, ErrStat, ErrMsg)
END SUBROUTINE SD_AB4

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
!! differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   Adams-Bashforth Predictor:
!!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!   Adams-Moulton Corrector:
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!!
!!  See, e.g.,
!!     - http://en.wikipedia.org/wiki/Linear_multistep_method
!!     - K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE SD_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other states at t on input at t + dt on output
      TYPE(SD_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ! local variables
      TYPE(SD_InputType)            :: u_interp        ! Continuous states at t
      TYPE(SD_ContinuousStateType)  :: x_pred          ! Continuous states at t
      TYPE(SD_ContinuousStateType)  :: xdot_pred       ! Continuous states at t

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL SD_CopyContState(x, x_pred, MESH_NEWCOPY, ErrStat, ErrMsg) !initialize x_pred      
      CALL SD_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, m, ErrStat, ErrMsg )

      if (n > 2) then
         CALL SD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg) ! make copy so that arrays/meshes get initialized/allocated for ExtrapInterp
         CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t + p%SDDeltaT, ErrStat, ErrMsg)

         CALL SD_CalcContStateDeriv(t + p%SDDeltaT, u_interp, p, x_pred, xd, z, OtherState, m, xdot_pred, ErrStat, ErrMsg ) ! initializes xdot_pred
         CALL SD_DestroyInput( u_interp, ErrStat, ErrMsg) ! local copy no longer needed

         x%qm    = x%qm    + (p%SDDeltaT / 24.) * ( 9. * xdot_pred%qm +  19. * OtherState%xdot(1)%qm - 5. * OtherState%xdot(2)%qm &
                                          + 1. * OtherState%xdot(3)%qm )
   
         x%qmdot = x%qmdot + (p%SDDeltaT / 24.) * ( 9. * xdot_pred%qmdot + 19. * OtherState%xdot(1)%qmdot - 5. * OtherState%xdot(2)%qmdot &
                                          + 1. * OtherState%xdot(3)%qmdot )
         CALL SD_DestroyContState( xdot_pred, ErrStat, ErrMsg) ! local copy no longer needed
      else
         x%qm    = x_pred%qm
         x%qmdot = x_pred%qmdot
      endif

      CALL SD_DestroyContState( x_pred, ErrStat, ErrMsg) ! local copy no longer needed
      
END SUBROUTINE SD_ABM4

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!!   Runge-Kutta." sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!!   Cambridge University Press, pp. 704-716, 1992.
SUBROUTINE SD_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other states at t on input at t + dt on output
      TYPE(SD_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ! local variables
      TYPE(SD_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(SD_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(SD_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(SD_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(SD_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(SD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(SD_InputType)                           :: u_interp    ! interpolated value of inputs 
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Initialize interim vars
      !bjj: the state type contains allocatable arrays, so we must first allocate space:
      CALL SD_CopyContState( x, k1,       MESH_NEWCOPY, ErrStat, ErrMsg )
      CALL SD_CopyContState( x, k2,       MESH_NEWCOPY, ErrStat, ErrMsg )
      CALL SD_CopyContState( x, k3,       MESH_NEWCOPY, ErrStat, ErrMsg )
      CALL SD_CopyContState( x, k4,       MESH_NEWCOPY, ErrStat, ErrMsg )
      CALL SD_CopyContState( x, x_tmp,    MESH_NEWCOPY, ErrStat, ErrMsg )
      
      ! interpolate u to find u_interp = u(t)
      CALL SD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg  )  ! we need to allocate input arrays/meshes before calling ExtrapInterp...     
      CALL SD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

      ! find xdot at t
      CALL SD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat, ErrMsg ) !initializes xdot
      k1%qm       = p%SDDeltaT * xdot%qm
      k1%qmdot    = p%SDDeltaT * xdot%qmdot
      x_tmp%qm    = x%qm    + 0.5 * k1%qm
      x_tmp%qmdot = x%qmdot + 0.5 * k1%qmdot
      ! interpolate u to find u_interp = u(t + dt/2)
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%SDDeltaT, ErrStat, ErrMsg)

      ! find xdot at t + dt/2
      CALL SD_CalcContStateDeriv( t + 0.5*p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )
      k2%qm    = p%SDDeltaT * xdot%qm
      k2%qmdot = p%SDDeltaT * xdot%qmdot
      x_tmp%qm    = x%qm    + 0.5 * k2%qm
      x_tmp%qmdot = x%qmdot + 0.5 * k2%qmdot

      ! find xdot at t + dt/2
      CALL SD_CalcContStateDeriv( t + 0.5*p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )
      k3%qm       = p%SDDeltaT * xdot%qm
      k3%qmdot    = p%SDDeltaT * xdot%qmdot
      x_tmp%qm    = x%qm    + k3%qm
      x_tmp%qmdot = x%qmdot + k3%qmdot
      ! interpolate u to find u_interp = u(t + dt)
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t + p%SDDeltaT, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL SD_CalcContStateDeriv( t + p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )
      k4%qm    = p%SDDeltaT * xdot%qm
      k4%qmdot = p%SDDeltaT * xdot%qmdot
      x%qm     = x%qm    +  ( k1%qm    + 2. * k2%qm    + 2. * k3%qm    + k4%qm    ) / 6.
      x%qmdot  = x%qmdot +  ( k1%qmdot + 2. * k2%qmdot + 2. * k3%qmdot + k4%qmdot ) / 6.

      CALL CleanUp()
      
CONTAINS       

   SUBROUTINE CleanUp()
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)
      CALL SD_DestroyContState( xdot,     ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )
      CALL SD_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
   END SUBROUTINE CleanUp            
      
END SUBROUTINE SD_RK4

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the 2nd-order Adams-Moulton Implicit Method (AM2,Trapezoidal rule) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 =  f(t       , x_t         )
!!        k2 =  f(t + dt  , x_t+dt      )
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) =x_n+1 = x_t + deltat/2*(k1 + k2) + O(dt^3)
!!   Now this can be re-written as: 0=Z(x_n+1) = x_n - x_n+1 +dt/2 *(f_n + f_n+1) = 0
!!         f_n= A*x_n + B*u_n + Fx  from Eq. 1.12 of the manual
!!         So to solve this linear system, I can just use x(k)=x(k-1) -J^-1 * Z(x(k-1))  (this is a simple root solver of the linear equation)
!!         with J=dZ/dx_n+1 = -I +dt/2*A 
!!
!!   Thus x_n+1 = x_n - J^-1 *dt/2 * (2*A*x_n + B *(u_n + u_n+1) +2*Fx)
!!  or    J*( x_n - x_n+1 ) = dt * ( A*x_n +  B *(u_n + u_n+1)/2 + Fx)
SUBROUTINE SD_AM2( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   REAL(DbKi),                     INTENT(IN   )   :: t              !< Current simulation time in seconds
   INTEGER(IntKi),                 INTENT(IN   )   :: n              !< time step number
   TYPE(SD_InputType),             INTENT(INOUT)   :: u(:)           !< Inputs at t
   REAL(DbKi),                     INTENT(IN   )   :: utimes(:)      !< times of input
   TYPE(SD_ParameterType),         INTENT(IN   )   :: p              !< Parameters
   TYPE(SD_ContinuousStateType),   INTENT(INOUT)   :: x              !< Continuous states at t on input at t + dt on output
   TYPE(SD_DiscreteStateType),     INTENT(IN   )   :: xd             !< Discrete states at t
   TYPE(SD_ConstraintStateType),   INTENT(IN   )   :: z              !< Constraint states at t (possibly a guess)
   TYPE(SD_OtherStateType),        INTENT(INOUT)   :: OtherState     !< Other states at t on input at t + dt on output
   TYPE(SD_MiscVarType),           INTENT(INOUT)   :: m              !< Misc/optimization variables
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   ! local variables
   TYPE(SD_InputType)                              :: u_interp       ! interpolated value of inputs 
   REAL(ReKi)                                      :: junk2(2*p%qml) !temporary states (qm and qmdot only)
   REAL(ReKi)                                      :: udot_TP2(6)    ! temporary copy of udot_TP
   REAL(ReKi)                                      :: udotdot_TP2(6) ! temporary copy of udotdot_TP
   REAL(ReKi)                                      :: UFL2(p%DOFL)   ! temporary copy of UFL
   INTEGER(IntKi)                                  :: ErrStat2
   CHARACTER(ErrMsgLen)                            :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = "" 

   ! Initialize interim vars
   CALL SD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2,ErrMsg2);CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
         
   !Start by getting u_n and u_n+1 
   ! interpolate u to find u_interp = u(t) = u_n     
   CALL SD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
   m%udotdot_TP = (/u_interp%TPMesh%TranslationAcc(:,1), u_interp%TPMesh%RotationAcc(:,1)/)
   m%udot_TP    = (/u_interp%TPMesh%TranslationVel(:,1), u_interp%TPMesh%RotationVel(:,1)/)
   CALL ConstructUFL( u_interp, p, m%UFL )     
                
   ! extrapolate u to find u_interp = u(t + dt)=u_n+1
   CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t+p%SDDeltaT, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
   udotdot_TP2 = (/u_interp%TPMesh%TranslationAcc(:,1), u_interp%TPMesh%RotationAcc(:,1)/)
   udot_TP2    = (/u_interp%TPMesh%TranslationVel(:,1), u_interp%TPMesh%RotationVel(:,1)/)
   CALL ConstructUFL( u_interp, p, UFL2 )     
   
   ! calculate (u_n + u_n+1)/2
   udot_TP2    = 0.5_ReKi * ( udot_TP2    + m%udot_TP )
   udotdot_TP2 = 0.5_ReKi * ( udotdot_TP2 + m%udotdot_TP )
   UFL2        = 0.5_ReKi * ( UFL2        + m%UFL        )
          
   ! set junk2 = dt * ( A*x_n +  B *(u_n + u_n+1)/2 + Fx)   
   junk2(      1:  p%qml)=p%SDDeltaT * x%qmdot                                                                                                   !upper portion of array
   junk2(1+p%qml:2*p%qml)=p%SDDeltaT * (p%NOmegaM2*x%qm + p%N2OmegaMJDamp*x%qmdot - matmul(p%CMM , x%qmdot)  &
                                                        - matmul(p%CMB, udot_TP2) - matmul(p%MMB, udotdot_TP2)  + matmul(UFL2,p%PhiM  ) + p%FX)  !lower portion of array
   ! note: matmul(UFL2,p%PhiM  ) = matmul(p%PhiM_T,UFL2) because UFL2 is 1-D

   IF (p%SeismicInp) THEN
      junk2(1+p%qml:2*p%qml) = junk2(1+p%qml:2*p%qml) + m%Ug * p%FMIMKP + m%Udotg * p%FMIMCP - m%Uddotg * p%FMIMMP + m%PHIg * p%FMIMKP_PHI + m%PHIdotg * p%FMIMCP_PHI - m%PHIddotg * p%FRIMMP_PHI + m%Vg * p%FMIMKP_V + m%Vdotg * p%FMIMCP_V - m%Vddotg * p%FMIMMP_V
   ENDIF
             
   !....................................................
   ! Solve for junk2: (equivalent to junk2= matmul(p%AM2InvJac,junk2)
   ! J*( x_n - x_n+1 ) = dt * ( A*x_n +  B *(u_n + u_n+1)/2 + Fx)
   !....................................................   
   CALL LAPACK_getrs( TRANS='N',N=SIZE(p%AM2Jac,1),A=p%AM2Jac,IPIV=p%AM2JacPiv, B=junk2, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
      !IF ( ErrStat >= AbortErrLev ) RETURN
      
   ! after the LAPACK solve, junk2 = ( x_n - x_n+1 ); so now we can solve for x_n+1:
   x%qm    = x%qm    - junk2(      1:  p%qml)
   x%qmdot = x%qmdot - junk2(p%qml+1:2*p%qml)
     
   ! clean up temporary variable(s)
   CALL SD_DestroyInput(  u_interp, ErrStat, ErrMsg )
   
END SUBROUTINE SD_AM2

!------------------------------------------------------------------------------------------------------
!> Perform Craig Bampton reduction
SUBROUTINE Craig_Bampton(Init, p, CBparams, ErrStat, ErrMsg)
   TYPE(SD_InitType),     INTENT(INOUT)      :: Init        ! Input data for initialization routine
   TYPE(SD_ParameterType),INTENT(INOUT)      :: p           ! Parameters
   TYPE(CB_MatArrays),    INTENT(INOUT)      :: CBparams    ! CB parameters that will be passed out for summary file use 
   INTEGER(IntKi),        INTENT(  OUT)      :: ErrStat     ! Error status of the operation
   CHARACTER(*),          INTENT(  OUT)      :: ErrMsg      ! Error message if ErrStat /= ErrID_None   
   ! local variables
   REAL(ReKi), ALLOCATABLE  :: MRR(:, :)
   REAL(ReKi), ALLOCATABLE  :: MLL(:, :)
   REAL(ReKi), ALLOCATABLE  :: MRL(:, :)
   REAL(ReKi), ALLOCATABLE  :: CRR(:, :)
   REAL(ReKi), ALLOCATABLE  :: CLL(:, :)
   REAL(ReKi), ALLOCATABLE  :: CRL(:, :)
   REAL(ReKi), ALLOCATABLE  :: KRR(:, :)
   REAL(ReKi), ALLOCATABLE  :: KLL(:, :)
   REAL(ReKi), ALLOCATABLE  :: KRL(:, :)
   REAL(ReKi), ALLOCATABLE  :: FGR(:)
   REAL(ReKi), ALLOCATABLE  :: FGL(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMKP(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMCP(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMMP(:)  !AÑADO FM
   REAL(ReKi), ALLOCATABLE  :: FMIMKP(:) 
   REAL(ReKi), ALLOCATABLE  :: FMIMCP(:)
   REAL(ReKi), ALLOCATABLE  :: FMIMMP(:)  !AÑADO FM
   !AÑADO GIROS
   REAL(ReKi), ALLOCATABLE  :: FRIMKP_PHI(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMCP_PHI(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMMP_PHI(:)  
   REAL(ReKi), ALLOCATABLE  :: FMIMKP_PHI(:) 
   REAL(ReKi), ALLOCATABLE  :: FMIMCP_PHI(:)
   REAL(ReKi), ALLOCATABLE  :: FMIMMP_PHI(:) 
   ! AÑADO VERTICAL
   REAL(ReKi), ALLOCATABLE  :: FRIMKP_V(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMCP_V(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMMP_V(:)  
   REAL(ReKi), ALLOCATABLE  :: FMIMKP_V(:) 
   REAL(ReKi), ALLOCATABLE  :: FMIMCP_V(:)
   REAL(ReKi), ALLOCATABLE  :: FMIMMP_V(:)  
   REAL(ReKi), ALLOCATABLE  :: MBBb(:, :)
   REAL(ReKi), ALLOCATABLE  :: MBMb(:, :) 
   REAL(ReKi), ALLOCATABLE  :: KBBb(:, :)
   REAL(ReKi), ALLOCATABLE  :: CBBb(:, :)
   REAL(ReKi), ALLOCATABLE  :: CBMb(:, :)
   REAL(ReKi), ALLOCATABLE  :: CMMb(:, :)
   REAL(ReKi), ALLOCATABLE  :: PhiRb(:, :)   
   REAL(ReKi), ALLOCATABLE  :: PhiRbase(:, :)   
   REAL(ReKi), ALLOCATABLE  :: FGRb(:) 
   REAL(ReKi), ALLOCATABLE  :: FRIMKPb(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMCPb(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMMPb(:) !AÑADO FM   
   REAL(ReKi), ALLOCATABLE  :: FMIMKPb(:)
   REAL(ReKi), ALLOCATABLE  :: FMIMCPb(:)
   REAL(ReKi), ALLOCATABLE  :: FMIMMPb(:)  !AÑADO FM
   !AÑADO GIROS
   REAL(ReKi), ALLOCATABLE  :: FRIMKP_PHIb(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMCP_PHIb(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMMP_PHIb(:)   
   REAL(ReKi), ALLOCATABLE  :: FMIMKP_PHIb(:)
   REAL(ReKi), ALLOCATABLE  :: FMIMCP_PHIb(:)
   REAL(ReKi), ALLOCATABLE  :: FMIMMP_PHIb(:)  
   ! AÑADO VERTICAL
   REAL(ReKi), ALLOCATABLE  :: FRIMKP_Vb(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMCP_Vb(:)
   REAL(ReKi), ALLOCATABLE  :: FRIMMP_Vb(:)   
   REAL(ReKi), ALLOCATABLE  :: FMIMKP_Vb(:)
   REAL(ReKi), ALLOCATABLE  :: FMIMCP_Vb(:)
   REAL(ReKi), ALLOCATABLE  :: FMIMMP_Vb(:)  
   REAL(ReKi)               :: JDamping1 ! temporary storage for first element of JDamping array 
   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(ErrMsgLen)     :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! number of nodes:
   p%NNodes_I  = Init%NInterf                         ! Number of interface nodes
   p%NNodes_L  = Init%NNode - p%NReact - p%NNodes_I   ! Number of Interior nodes =(TDOF-DOFC-DOFI)/6 =  (6*Init%NNode - (p%NReact+p%NNodes_I)*6 ) / 6 = Init%NNode - p%NReact -p%NNodes_I

   !DOFS of interface
   !BJJ: TODO:  are these 6's actually TPdofL?   
   p%DOFI = p%NNodes_I*6
   p%DOFC = p%NReact*6
   p%DOFR = (p%NReact+p%NNodes_I)*6 ! = p%DOFC + p%DOFI
   p%DOFL = p%NNodes_L*6            ! = Init%TDOF - p%DOFR
   
            
   IF(Init%CBMod) THEN ! C-B reduction         
      ! check number of internal modes
      IF(p%Nmodes > p%DOFL) THEN
         CALL SetErrStat(ErrID_Fatal,'Number of internal modes is larger than maximum. ',ErrStat,ErrMsg,'Craig_Bampton')
         CALL CleanupCB()
         RETURN
      ENDIF
      
   ELSE ! full FEM 
      p%Nmodes = p%DOFL
      !Jdampings  need to be reallocated here because DOFL not known during Init
      !So assign value to one temporary variable
      JDamping1=Init%Jdampings(1)
      DEALLOCATE(Init%JDampings)
      CALL AllocAry( Init%JDampings, p%DOFL, 'Init%JDampings',  ErrStat2, ErrMsg2 ) ; if(Failed()) return
      Init%JDampings = JDamping1 ! set default values for all modes
      
   ENDIF   
   
   CBparams%DOFM = p%Nmodes  ! retained modes (all if no C-B reduction)
      
   ! matrix dimension paramters
   p%qmL    = p%Nmodes                       ! Length of 1/2 x array, x1 that is (note, do this after check if CBMod is true [Nmodes modified if CMBod is false])
   p%URbarL = p%DOFI !=p%NNodes_I*6          ! Length of URbar array, subarray of Y2  : THIS MAY CHANGE IF SOME DOFS ARE NOT CONSTRAINED       
   
      
   CALL AllocParameters(p, CBparams%DOFM, ErrStat2, ErrMsg2);                                  ; if (Failed()) return

   CALL AllocAry( MRR,             p%DOFR, p%DOFR,        'matrix MRR',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( MLL,             p%DOFL, p%DOFL,        'matrix MLL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( MRL,             p%DOFR, p%DOFL,        'matrix MRL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( CRR,             p%DOFR, p%DOFR,        'matrix CRR',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')   
   CALL AllocAry( CLL,             p%DOFL, p%DOFL,        'matrix CLL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( CRL,             p%DOFR, p%DOFL,        'matrix CRL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( KRR,             p%DOFR, p%DOFR,        'matrix KRR',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( KLL,             p%DOFL, p%DOFL,        'matrix KLL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( KRL,             p%DOFR, p%DOFL,        'matrix KRL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FGL,             p%DOFL,                'array FGL',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FGR,             p%DOFR,                'array FGR',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FMIMKP,          CBparams%DOFM,         'array FMIMKP,',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FRIMKP,          p%DOFI,                'array FRIMKP',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( FMIMCP,          CBparams%DOFM,         'array FMIMCP,',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FRIMCP,          p%DOFI,                'array FRIMCP',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( FMIMMP,          CBparams%DOFM,         'array FMIMMP,',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') !AÑADO FM  
   CALL AllocAry( FRIMMP,          p%DOFI,                'array FRIMMP',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') !AÑADO FM 
   ! AÑADO GIROS  
   CALL AllocAry( FMIMKP_PHI,          CBparams%DOFM,         'array FMIMKP_PHI,',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FRIMKP_PHI,          p%DOFI,                'array FRIMKP_PHI',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( FMIMCP_PHI,          CBparams%DOFM,         'array FMIMCP_PHI,',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FRIMCP_PHI,          p%DOFI,                'array FRIMCP_PHI',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( FMIMMP_PHI,          CBparams%DOFM,         'array FMIMMP_PHI,',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( FRIMMP_PHI,          p%DOFI,                'array FRIMMP_PHI',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   ! AÑADO VERTICAL
   CALL AllocAry( FMIMKP_V,          CBparams%DOFM,         'array FMIMKP_V,',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FRIMKP_V,          p%DOFI,                'array FRIMKP_V',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( FMIMCP_V,          CBparams%DOFM,         'array FMIMCP_V,',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FRIMCP_V,          p%DOFI,                'array FRIMCP_V',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( FMIMMP_V,          CBparams%DOFM,         'array FMIMMP_V,',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
   CALL AllocAry( FRIMMP_V,          p%DOFI,                'array FRIMMP_V',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')   
   CALL AllocAry( CBparams%MBB,    p%DOFR, p%DOFR,       'CBparams%MBB',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%MBM,    p%DOFR, CBparams%DOFM,'CBparams%MBM',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%KBB,    p%DOFR, p%DOFR,       'CBparams%KBB',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%CBB,    p%DOFR, p%DOFR,       'CBparams%CBB',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%CBM,    p%DOFR, CBparams%DOFM,'CBparams%CBM',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%CMM,    CBparams%DOFM, CBparams%DOFM,'CBparams%CMM',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')

   CALL AllocAry( CBparams%PhiL,   p%DOFL, p%DOFL,       'CBparams%PhiL',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%PhiR,   p%DOFL, p%DOFR,       'CBparams%PhiR',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%OmegaL, p%DOFL,               'CBparams%OmegaL', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%TI2,    p%DOFR, 6,            'CBparams%TI2',    ErrStat2, ErrMsg2 ); if(Failed()) return
   
   ! Set the index arrays p%IDI, p%IDR, p%IDL, p%IDC, and p%IDY. 
   CALL SetIndexArrays(Init, p, ErrStat2, ErrMsg2) ; if(Failed()) return

   ! Set MRR, MLL, MRL, CRR, CLL, CRL, KRR, KLL, KRL, FGR, FGL, based on
   !     Init%M, Init%K, and Init%FG data and indices p%IDR and p%IDL:
   CALL BreakSysMtrx(Init, p, MRR, MLL, MRL, CRR, CLL, CRL, KRR, KLL, KRL, FGR, FGL)   
      
   ! Set p%TI and CBparams%TI2
   CALL TrnsfTI(Init, p%TI, p%DOFI, p%IDI, CBparams%TI2, p%DOFR, p%IDR, ErrStat2, ErrMsg2); if(Failed()) return

   !................................
   ! Sets the following values, as documented in the SubDyn Theory Guide:
   !    CBparams%OmegaL (omega) and CBparams%PhiL from Eq. 2
   !    p%PhiL_T and p%PhiLInvOmgL2 for static improvement 
   !    CBparams%PhiR from Eq. 3
   !    CBparams%MBB, CBparams%MBM, and CBparams%KBB from Eq. 4.
   !................................
   CALL CBMatrix(MRR, MLL, MRL, CRR, CLL, CRL, KRR, KLL, KRL, CBparams%DOFM, Init, &  ! < inputs
                 CBparams%MBB, CBparams%MBM, CBparams%KBB,  CBparams%CBB, CBparams%CBM, CBparams%CMM, &  ! <- outputs
                 FRIMKP, FMIMKP, FRIMCP, FMIMCP, FRIMMP, FMIMMP, & !AÑADO FM
                 FRIMKP_PHI, FMIMKP_PHI, FRIMCP_PHI, FMIMCP_PHI, FRIMMP_PHI, FMIMMP_PHI, & !AÑADO GIROS
                 FRIMKP_V, FMIMKP_V, FRIMCP_V, FMIMCP_V, FRIMMP_V, FMIMMP_V, & ! AÑADO VERTICAL
                 CBparams%PhiL, CBparams%PhiR, CBparams%OmegaL, ErrStat2, ErrMsg2, p)  ! <- outputs (p is also input )
   if(Failed()) return
      
   ! to use a little less space, let's deallocate these arrays that we don't need anymore, then allocate the next set of temporary arrays:     
   IF(ALLOCATED(MRR)  ) DEALLOCATE(MRR) 
   IF(ALLOCATED(MLL)  ) DEALLOCATE(MLL) 
   IF(ALLOCATED(MRL)  ) DEALLOCATE(MRL) 
   IF(ALLOCATED(CRR)  ) DEALLOCATE(CRR) 
   IF(ALLOCATED(CLL)  ) DEALLOCATE(CLL) 
   IF(ALLOCATED(CRL)  ) DEALLOCATE(CRL) 
   IF(ALLOCATED(KRR)  ) DEALLOCATE(KRR) 
   IF(ALLOCATED(KLL)  ) DEALLOCATE(KLL) 
   IF(ALLOCATED(KRL)  ) DEALLOCATE(KRL) 

   ! "b" stands for "bar"; "t" stands for "tilde"
   CALL AllocAry( MBBb,  p%DOFI, p%DOFI,               'matrix MBBb',  ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( MBmb,  p%DOFI, CBparams%DOFM,        'matrix MBmb',  ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( KBBb,  p%DOFI, p%DOFI,               'matrix KBBb',  ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( CBBb,  p%DOFI, p%DOFI,               'matrix CBBb',  ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( CBmb,  p%DOFI, CBparams%DOFM,        'matrix CBmb',  ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( CMMb,  CBparams%DOFM, CBparams%DOFM, 'matrix CMMb',  ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( PhiRb, p%DOFL, p%DOFI,               'matrix PhiRb', ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( PhiRbase, p%DOFL, p%DOFR-p%DOFI,     'matrix PhiRbase', ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( FGRb,  p%DOFI,                       'array FGRb',   ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FRIMKPb,p%DOFI,                       'array FRIMKPb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FRIMCPb,p%DOFI,                       'array FRIMCPb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FRIMMPb,p%DOFI,                       'array FRIMMPb',ErrStat2, ErrMsg2 ); if (Failed()) return  !AÑADO FM 
   CALL AllocAry(FMIMKPb,CBparams%DOFM,                'array FMIMKPb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FMIMCPb,CBparams%DOFM,                'array FMIMCPb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FMIMMPb,CBparams%DOFM,                'array FMIMMPb',ErrStat2, ErrMsg2 ); if (Failed()) return  !AÑADO FM  
   !AÑADO GIROS
   CALL AllocAry(FRIMKP_PHIb,p%DOFI,                       'array FRIMKP_PHIb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FRIMCP_PHIb,p%DOFI,                       'array FRIMCP_PHIb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FRIMMP_PHIb,p%DOFI,                       'array FRIMMP_PHIb',ErrStat2, ErrMsg2 ); if (Failed()) return  
   CALL AllocAry(FMIMKP_PHIb,CBparams%DOFM,                'array FMIMKP_PHIb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FMIMCP_PHIb,CBparams%DOFM,                'array FMIMCP_PHIb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FMIMMP_PHIb,CBparams%DOFM,                'array FMIMMP_PHIb',ErrStat2, ErrMsg2 ); if (Failed()) return  
   !AÑADO VERTICAL
   CALL AllocAry(FRIMKP_Vb,p%DOFI,                       'array FRIMKP_Vb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FRIMCP_Vb,p%DOFI,                       'array FRIMCP_Vb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FRIMMP_Vb,p%DOFI,                       'array FRIMMP_Vb',ErrStat2, ErrMsg2 ); if (Failed()) return  
   CALL AllocAry(FMIMKP_Vb,CBparams%DOFM,                'array FMIMKP_Vb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FMIMCP_Vb,CBparams%DOFM,                'array FMIMCP_Vb',ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry(FMIMMP_Vb,CBparams%DOFM,                'array FMIMMP_Vb',ErrStat2, ErrMsg2 ); if (Failed()) return    
   
   !................................
   ! Convert CBparams%MBB , CBparams%MBM , CBparams%KBB , CBparams%CBB , CBparams%CBM , CBparams%CMM , CBparams%PhiR , FGR to
   !                  MBBb,          MBMb,          KBBb,          CBBb,          CBMb,          CMMb,          PHiRb, FGRb
   ! (throw out rows/columns of first matrices to create second matrices)
   !................................
   CALL CBApplyConstr(p%DOFI, p%DOFR, CBparams%DOFM,  p%DOFL,  &
                      CBparams%MBB , CBparams%MBM , CBparams%KBB ,                             &
                      CBparams%CBB , CBparams%CBM , CBparams%CMM , CBparams%PhiR , FGR ,       &
                            FRIMKP ,       FRIMCP ,       FRIMMP ,                             & !AÑADO FRM
                            FMIMKP ,       FMIMCP ,       FMIMMP ,                             & !AÑADO FM
                            FRIMKPb,       FRIMCPb,       FRIMMPb,                             & !AÑADO FMb
                            FMIMKPb,       FMIMCPb,       FMIMMPb,                             & !AÑADO FMb
                            !AÑADO GIROS
                            FRIMKP_PHI ,       FRIMCP_PHI ,       FRIMMP_PHI ,                 & 
                            FMIMKP_PHI ,       FMIMCP_PHI ,       FMIMMP_PHI ,                 & 
                            FRIMKP_PHIb,       FRIMCP_PHIb,       FRIMMP_PHIb,                 & 
                            FMIMKP_PHIb,       FMIMCP_PHIb,       FMIMMP_PHIb,                 & 
                            !AÑADO VERTICAL
                            FRIMKP_V ,       FRIMCP_V ,       FRIMMP_V ,                       & 
                            FMIMKP_V ,       FMIMCP_V ,       FMIMMP_V ,                       & 
                            FRIMKP_Vb,       FRIMCP_Vb,       FRIMMP_Vb,                       & 
                            FMIMKP_Vb,       FMIMCP_Vb,       FMIMMP_Vb,                       &                             
                               MBBb,          MBMb,          KBBb,                             &                             
                               CBBb,          CBMb,          CMMb,          PHiRb, PhiRbase, FGRb, p)

   !................................
   ! set values needed to calculate outputs and update states:
   !................................
   CALL SetParameters(Init, p, MBBb, MBmb, KBBb, CBBb, CBmb, CMMb, FGRb, &
                      FRIMKPb, FRIMCPb,FRIMMPb, FMIMKPb, FMIMCPb, FMIMMPb,FRIMKP_PHIb, FRIMCP_PHIb,FRIMMP_PHIb, FMIMKP_PHIb, FMIMCP_PHIb, FMIMMP_PHIb, &
                      FRIMKP_Vb, FRIMCP_Vb,FRIMMP_Vb, FMIMKP_Vb, FMIMCP_Vb, FMIMMP_Vb, PhiRb, PhiRbase, CBparams%OmegaL, FGL, CBparams%PhiL, ErrStat2, ErrMsg2)  !AÑADO FM,  
                      !AÑADIR GIRO BASE  !AÑADIR VERTICAL BASE                 
   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'Craig_Bampton')
      
   CALL CleanUpCB()

contains

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUpCB()
   end function Failed

   subroutine CleanUpCB()
      IF(ALLOCATED(MRR)  ) DEALLOCATE(MRR) 
      IF(ALLOCATED(MLL)  ) DEALLOCATE(MLL) 
      IF(ALLOCATED(MRL)  ) DEALLOCATE(MRL) 
      IF(ALLOCATED(CRR)  ) DEALLOCATE(CRR) 
      IF(ALLOCATED(CLL)  ) DEALLOCATE(CLL) 
      IF(ALLOCATED(CRL)  ) DEALLOCATE(CRL)
      IF(ALLOCATED(KRR)  ) DEALLOCATE(KRR) 
      IF(ALLOCATED(KLL)  ) DEALLOCATE(KLL) 
      IF(ALLOCATED(KRL)  ) DEALLOCATE(KRL) 
      IF(ALLOCATED(FGL)  ) DEALLOCATE(FGL) 
      IF(ALLOCATED(FGR)  ) DEALLOCATE(FGR) 
      IF(ALLOCATED(FRIMKP)) DEALLOCATE(FRIMKP)
      IF(ALLOCATED(FMIMKP)) DEALLOCATE(FMIMKP)
      IF(ALLOCATED(FRIMCP)) DEALLOCATE(FRIMCP)
      IF(ALLOCATED(FMIMCP)) DEALLOCATE(FMIMCP)
      IF(ALLOCATED(FRIMMP)) DEALLOCATE(FRIMMP) !AÑADO FM
      IF(ALLOCATED(FMIMMP)) DEALLOCATE(FMIMMP) !AÑADO FM     
      !AÑADO GIRO
      IF(ALLOCATED(FRIMKP_PHI)) DEALLOCATE(FRIMKP_PHI)
      IF(ALLOCATED(FMIMKP_PHI)) DEALLOCATE(FMIMKP_PHI)
      IF(ALLOCATED(FRIMCP_PHI)) DEALLOCATE(FRIMCP_PHI)
      IF(ALLOCATED(FMIMCP_PHI)) DEALLOCATE(FMIMCP_PHI)
      IF(ALLOCATED(FRIMMP_PHI)) DEALLOCATE(FRIMMP_PHI) 
      IF(ALLOCATED(FMIMMP_PHI)) DEALLOCATE(FMIMMP_PHI)    
      ! AÑADO VERTICAL
      IF(ALLOCATED(FRIMKP_V)) DEALLOCATE(FRIMKP_V)
      IF(ALLOCATED(FMIMKP_V)) DEALLOCATE(FMIMKP_V)
      IF(ALLOCATED(FRIMCP_V)) DEALLOCATE(FRIMCP_V)
      IF(ALLOCATED(FMIMCP_V)) DEALLOCATE(FMIMCP_V)
      IF(ALLOCATED(FRIMMP_V)) DEALLOCATE(FRIMMP_V) 
      IF(ALLOCATED(FMIMMP_V)) DEALLOCATE(FMIMMP_V)       
      IF(ALLOCATED(MBBb) ) DEALLOCATE(MBBb) 
      IF(ALLOCATED(MBmb) ) DEALLOCATE(MBmb) 
      IF(ALLOCATED(KBBb) ) DEALLOCATE(KBBb) 
      IF(ALLOCATED(CBBb) ) DEALLOCATE(CBBb) 
      IF(ALLOCATED(CBmb) ) DEALLOCATE(CBmb) 
      IF(ALLOCATED(CMMb) ) DEALLOCATE(CMMb) 
      IF(ALLOCATED(PhiRb)) DEALLOCATE(PhiRb) 
      IF(ALLOCATED(PhiRbase)) DEALLOCATE(PhiRbase) 
      IF(ALLOCATED(FGRb) ) DEALLOCATE(FGRb) 
      IF(ALLOCATED(FRIMKPb)) DEALLOCATE(FRIMKPb)
      IF(ALLOCATED(FMIMKPb)) DEALLOCATE(FMIMKPb)
      IF(ALLOCATED(FRIMCPb)) DEALLOCATE(FRIMCPb)
      IF(ALLOCATED(FMIMCPb)) DEALLOCATE(FMIMCPb)  
      IF(ALLOCATED(FRIMMPb)) DEALLOCATE(FRIMMPb) !AÑADO FM
      IF(ALLOCATED(FMIMMPb)) DEALLOCATE(FMIMMPb) !AÑADO FM 
      !AÑADO GIROb
      IF(ALLOCATED(FRIMKP_PHIb)) DEALLOCATE(FRIMKP_PHIb)
      IF(ALLOCATED(FMIMKP_PHIb)) DEALLOCATE(FMIMKP_PHIb)
      IF(ALLOCATED(FRIMCP_PHIb)) DEALLOCATE(FRIMCP_PHIb)
      IF(ALLOCATED(FMIMCP_PHIb)) DEALLOCATE(FMIMCP_PHIb)  
      IF(ALLOCATED(FRIMMP_PHIb)) DEALLOCATE(FRIMMP_PHIb) 
      IF(ALLOCATED(FMIMMP_PHIb)) DEALLOCATE(FMIMMP_PHIb) 
      !AÑADO VERTICALb
      IF(ALLOCATED(FRIMKP_Vb)) DEALLOCATE(FRIMKP_Vb)
      IF(ALLOCATED(FMIMKP_Vb)) DEALLOCATE(FMIMKP_Vb)
      IF(ALLOCATED(FRIMCP_Vb)) DEALLOCATE(FRIMCP_Vb)
      IF(ALLOCATED(FMIMCP_Vb)) DEALLOCATE(FMIMCP_Vb)  
      IF(ALLOCATED(FRIMMP_Vb)) DEALLOCATE(FRIMMP_Vb) 
      IF(ALLOCATED(FMIMMP_Vb)) DEALLOCATE(FMIMMP_Vb)   
                
   end subroutine CleanUpCB

END SUBROUTINE Craig_Bampton 

!------------------------------------------------------------------------------------------------------
!>
SUBROUTINE BreakSysMtrx(Init, p, MRR, MLL, MRL, CRR, CLL, CRL, KRR, KLL, KRL, FGR, FGL   )
   TYPE(SD_InitType),      INTENT(IN   )  :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType), INTENT(IN   )  :: p  
   REAL(ReKi),             INTENT(  OUT)  :: MRR(p%DOFR, p%DOFR)
   REAL(ReKi),             INTENT(  OUT)  :: MLL(p%DOFL, p%DOFL) 
   REAL(ReKi),             INTENT(  OUT)  :: MRL(p%DOFR, p%DOFL)
   REAL(ReKi),             INTENT(  OUT)  :: CRR(p%DOFR, p%DOFR)
   REAL(ReKi),             INTENT(  OUT)  :: CLL(p%DOFL, p%DOFL) 
   REAL(ReKi),             INTENT(  OUT)  :: CRL(p%DOFR, p%DOFL)
   REAL(ReKi),             INTENT(  OUT)  :: KRR(p%DOFR, p%DOFR)
   REAL(ReKi),             INTENT(  OUT)  :: KLL(p%DOFL, p%DOFL)
   REAL(ReKi),             INTENT(  OUT)  :: KRL(p%DOFR, p%DOFL)
   REAL(ReKi),             INTENT(  OUT)  :: FGR(p%DOFR)
   REAL(ReKi),             INTENT(  OUT)  :: FGL(p%DOFL)
   ! local variables
   INTEGER(IntKi)          :: I, J, II, JJ
   
   DO I = 1, p%DOFR   !Boundary DOFs
      II = p%IDR(I)
      FGR(I) = Init%FG(II)
      DO J = 1, p%DOFR
         JJ = p%IDR(J)
         MRR(I, J) = Init%M(II, JJ)
         CRR(I, J) = Init%C(II, JJ)
         KRR(I, J) = Init%K(II, JJ)
      ENDDO
   ENDDO
   
   DO I = 1, p%DOFL
      II = p%IDL(I)
      FGL(I) = Init%FG(II)
      DO J = 1, p%DOFL
         JJ = p%IDL(J)
         MLL(I, J) = Init%M(II, JJ)
         CLL(I, J) = Init%C(II, JJ)
         KLL(I, J) = Init%K(II, JJ)
      ENDDO
   ENDDO
   
   DO I = 1, p%DOFR
      II = p%IDR(I)
      DO J = 1, p%DOFL
         JJ = p%IDL(J)
         MRL(I, J) = Init%M(II, JJ)
         CRL(I, J) = Init%C(II, JJ)
         KRL(I, J) = Init%K(II, JJ)   !Note KRL and MRL are getting data from a constraint-applied formatted M and K (i.e. Mbar and Kbar) this may not be legit!! RRD
      ENDDO                           !I think this is fixed now since the constraint application occurs later
   ENDDO
      
END SUBROUTINE BreakSysMtrx

!------------------------------------------------------------------------------------------------------
!> Sets the CB values, as documented in the SubDyn Theory Guide:
! OmegaL (omega) and PhiL from Eq. 2
! p%PhiL_T and p%PhiLInvOmgL2 for static improvement (will be added to theory guide later?)
! PhiR from Eq. 3
! MBB, MBM, and KBB from Eq. 4.
!................................
SUBROUTINE CBMatrix( MRR, MLL, MRL, CRR, CLL, CRL, KRR, KLL, KRL, DOFM, Init, & !AÑADIR MIB, MLB 
                     MBB, MBM, KBB, CBB, CBM, CMM, FRIMKP, FMIMKP, FRIMCP, FMIMCP, FRIMMP, FMIMMP, &
                     FRIMKP_PHI, FMIMKP_PHI, FRIMCP_PHI, FMIMCP_PHI, FRIMMP_PHI, FMIMMP_PHI,       &
                     FRIMKP_V, FMIMKP_V, FRIMCP_V, FMIMCP_V, FRIMMP_V, FMIMMP_V,                   &
                      PhiL, PhiR, OmegaL, ErrStat, ErrMsg, p) ! AÑADIR GIRO, AÑADIR VERTICAL

   TYPE(SD_InitType),      INTENT(IN)    :: Init
   TYPE(SD_ParameterType), INTENT(INOUT) :: p  
   INTEGER(IntKi),         INTENT(  in)  :: DOFM
   REAL(ReKi),             INTENT(  IN)  :: MRR( p%DOFR, p%DOFR)
   REAL(ReKi),             INTENT(  IN)  :: MLL( p%DOFL, p%DOFL) 
   REAL(ReKi),             INTENT(  IN)  :: MRL( p%DOFR, p%DOFL)
   REAL(ReKi),             INTENT(  IN)  :: CRR( p%DOFR, p%DOFR)
   REAL(ReKi),             INTENT(  IN)  :: CLL( p%DOFL, p%DOFL) 
   REAL(ReKi),             INTENT(  IN)  :: CRL( p%DOFR, p%DOFL)
   REAL(ReKi),             INTENT(  IN)  :: KRR( p%DOFR, p%DOFR)
   REAL(ReKi),             INTENT(INOUT) :: KLL( p%DOFL, p%DOFL)  ! on exit, it has been factored (otherwise not changed)
   REAL(ReKi),             INTENT(  IN)  :: KRL( p%DOFR, p%DOFL)
   REAL(ReKi),             INTENT(INOUT) :: MBB( p%DOFR, p%DOFR)
   REAL(ReKi),             INTENT(INOUT) :: MBM( p%DOFR,   DOFM)
   REAL(ReKi),             INTENT(INOUT) :: KBB( p%DOFR, p%DOFR)
   REAL(ReKi),             INTENT(INOUT) :: CBB( p%DOFR, p%DOFR)
   REAL(ReKi),             INTENT(INOUT) :: CBM( p%DOFR,   DOFM)
   REAL(ReKi),             INTENT(INOUT) :: CMM(   DOFM,   DOFM)
  ! REAL(ReKi),             INTENT(INOUT) :: MIb( p%DOFI,P%DOFR - P%DOFI)
  ! REAL(ReKi),             INTENT(INOUT) :: MMb(   DOFM,P%DOFR - P%DOFI)     
   REAL(ReKi),             INTENT(INOUT) :: FRIMKP (   p%DOFI)
   REAL(ReKi),             INTENT(INOUT) :: FRIMCP (   p%DOFI)
   REAL(ReKi),             INTENT(INOUT) :: FRIMMP (   p%DOFI) !AÑADO FM  
   REAL(ReKi),             INTENT(INOUT) :: FMIMKP (   DOFM  )
   REAL(ReKi),             INTENT(INOUT) :: FMIMCP (   DOFM  )
   REAL(ReKi),             INTENT(INOUT) :: FMIMMP (   DOFM  ) !AÑADO FM
   ! AÑADO GIROS
   REAL(ReKi),             INTENT(INOUT) :: FRIMKP_PHI (   p%DOFI)
   REAL(ReKi),             INTENT(INOUT) :: FRIMCP_PHI (   p%DOFI)
   REAL(ReKi),             INTENT(INOUT) :: FRIMMP_PHI (   p%DOFI) 
   REAL(ReKi),             INTENT(INOUT) :: FMIMKP_PHI (   DOFM  )
   REAL(ReKi),             INTENT(INOUT) :: FMIMCP_PHI (   DOFM  )
   REAL(ReKi),             INTENT(INOUT) :: FMIMMP_PHI (   DOFM  )    
   ! VERTICAL
   REAL(ReKi),             INTENT(INOUT) :: FRIMKP_V (   p%DOFI)
   REAL(ReKi),             INTENT(INOUT) :: FRIMCP_V (   p%DOFI)
   REAL(ReKi),             INTENT(INOUT) :: FRIMMP_V (   p%DOFI) 
   REAL(ReKi),             INTENT(INOUT) :: FMIMKP_V (   DOFM  )
   REAL(ReKi),             INTENT(INOUT) :: FMIMCP_V (   DOFM  )
   REAL(ReKi),             INTENT(INOUT) :: FMIMMP_V (   DOFM  ) 
   REAL(ReKi),             INTENT(INOUT) :: PhiR(p%DOFL, p%DOFR)   
   REAL(ReKi),             INTENT(INOUT) :: PhiL(p%DOFL, p%DOFL)    !used to be PhiM(DOFL,DOFM), now it is more generic
   REAL(ReKi),             INTENT(INOUT) :: OmegaL(p%DOFL)   !used to be omegaM only   ! Eigenvalues
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! LOCAL VARIABLES
   REAL(ReKi) , allocatable               :: Mu(:, :)          ! matrix for normalization Mu(p%DOFL, p%DOFL) [bjj: made allocatable to try to avoid stack issues]
   REAL(ReKi) , allocatable               :: Temp(:, :)        ! temp matrix for intermediate steps [bjj: made allocatable to try to avoid stack issues]
   REAL(ReKi) , allocatable               :: OmegaDamp2(:, :)
   REAL(ReKi) , allocatable               :: NOmegaM2(:)
   REAL(ReKi) , allocatable               :: InvPhiM(:, :)
   REAL(ReKi) , allocatable               :: PhiR_T_MLL(:,:)   ! PhiR_T_MLL(p%DOFR,p%DOFL) = transpose of PhiR * MLL (temporary storage)
   REAL(ReKi) , allocatable               :: PhiR_T_CLL(:,:)   ! PhiR_T_CLL(p%DOFR,p%DOFL) = transpose of PhiR * CLL (temporary storage)
   REAL(ReKi) , allocatable               :: RR(:)             ! RR(p%DOFR) , influece vector on R dofs
   REAL(ReKi) , allocatable               :: RI(:)             ! RI(p%DOFR) , influece vector on I dofs
   REAL(ReKi) , allocatable               :: RL(:)             ! RR(p%DOFL) , influece vector on L dofs
   REAL(ReKi) , allocatable               :: RB(:)             ! RR(p%DOFL) , influece vector on b dofs   
   REAL(ReKi) , allocatable               :: RR_PHI(:)             ! RR(p%DOFR) , influece vector on R dofs
   REAL(ReKi) , allocatable               :: RI_PHI(:)             ! RR(p%DOFR) , influece vector on I dofs
   REAL(ReKi) , allocatable               :: RL_PHI(:)             ! RR(p%DOFL) , influece vector on L dofs
   REAL(ReKi) , allocatable               :: RB_PHI(:)             ! RR(p%DOFL) , influece vector on b dofs  
   REAL(ReKi) , allocatable               :: RR_V(:)             ! RR(p%DOFR) , influece vector on R dofs
   REAL(ReKi) , allocatable               :: RI_V(:)             ! RR(p%DOFR) , influece vector on I dofs
   REAL(ReKi) , allocatable               :: RL_V(:)             ! RR(p%DOFL) , influece vector on L dofs
   REAL(ReKi) , allocatable               :: RB_V(:)             ! RR(p%DOFL) , influece vector on b dofs     
   REAL(ReKi) , allocatable               :: KLLcopy(:,:)     ! we need a copy to use after KLL has been factored
   INTEGER                                :: I, J !, lwork !counter, and varibales for inversion routines
   INTEGER                                :: DOFvar !placeholder used to get both PhiL or PhiM into 1 process
   INTEGER                                :: ipiv(p%DOFL) !the integer vector ipvt of length min(m,n), containing the pivot indices. 
                                                       !Returned as: a one-dimensional array of (at least) length min(m,n), containing integers,
                                                       !where 1 <= less than or equal to ipvt(i) <= less than or equal to m.
   INTEGER(IntKi)                         :: ErrStat2                                                                    
   CHARACTER(ErrMsgLen)                   :: ErrMsg2
   CHARACTER(*), PARAMETER                :: RoutineName = 'CBMatrix'
                                                       
   ErrStat = ErrID_None 
   ErrMsg  = ''
   
   CALL WrScr('   Calculating Internal Modal Eigenvectors')
        
   IF (p%SttcSolve) THEN ! STATIC TREATMENT IMPROVEMENT
      DOFvar=p%DOFL
   ELSE
      DOFvar=DOFM !Initialize for normal cases, dynamic only      
   ENDIF  

   CALL AllocAry( KLLcopy , p%DOFL , p%DOFL , 'KLLcopy' , ErrStat2 , ErrMsg2); if(Failed()) return
   KLLcopy = KLL     

   !....................................................
   ! Set OmegaL and PhiL from Eq. 2
   !....................................................
   IF ( DOFvar > 0 ) THEN ! Only time this wouldn't happen is if no modes retained and no static improvement...
      CALL EigenSolve(KLL, MLL, p%DOFL, DOFvar, .False.,Init,p, PhiL(:,1:DOFvar), OmegaL(1:DOFvar),  ErrStat2, ErrMsg2); if(Failed()) return

      ! --- Normalize PhiL
      ! bjj: break up this equation to avoid as many tenporary variables on the stack
      ! MU = MATMUL ( MATMUL( TRANSPOSE(PhiL), MLL ), PhiL )
      CALL AllocAry( Temp , p%DOFL , p%DOFL , 'Temp' , ErrStat2 , ErrMsg2); if(Failed()) return
      CALL AllocAry( MU   , p%DOFL , p%DOFL , 'Mu'   , ErrStat2 , ErrMsg2); if(Failed()) return
      MU   = TRANSPOSE(PhiL)
      Temp = MATMUL( MU, MLL )
      MU   = MATMUL( Temp, PhiL )
      DEALLOCATE(Temp)
      ! PhiL = MATMUL( PhiL, MU2 )  !this is the nondimensionalization (MU2 is diagonal)   
      DO I = 1, DOFvar
         PhiL(:,I) = PhiL(:,I) / SQRT( MU(I, I) )
      ENDDO    
      DO I=DOFvar+1, p%DOFL !loop done only if .not. p%SttcSolve .and. DOFM < p%DOFL (and actually, in that case, these values aren't used anywhere anyway)
         PhiL(:,I) = 0.0_ReKi
         OmegaL(I) = 0.0_ReKi
      END DO     
      DEALLOCATE(MU)
      
      !....................................................
      ! Set p%PhiL_T and p%PhiLInvOmgL2 for static improvement
      !....................................................
      IF (p%SttcSolve) THEN   
         p%PhiL_T=TRANSPOSE(PhiL) !transpose of PhiL for static improvement
         DO I = 1, p%DOFL
            p%PhiLInvOmgL2(:,I) = PhiL(:,I)* (1./OmegaL(I)**2)
         ENDDO 
      END IF
      
   ! ELSE .not. p%SttcSolve .and. DOFM < p%DOFL (in this case, PhiL, OmegaL aren't used)      
   END IF
      
      
   !....................................................
   ! Set PhiR from Eq. 3:
   !....................................................   
   ! now factor KLL to compute PhiR: KLL*PhiR=-TRANSPOSE(KRL)
   ! ** note this must be done after EigenSolve() because it modifies KLL **
   CALL LAPACK_getrf( p%DOFL, p%DOFL, KLL, ipiv, ErrStat2, ErrMsg2); if(Failed()) return
   
   PhiR = -1.0_ReKi * TRANSPOSE(KRL) !set "b" in Ax=b  (solve KLL * PhiR = - TRANSPOSE( KRL ) for PhiR)
   CALL LAPACK_getrs( TRANS='N',N=p%DOFL,A=KLL,IPIV=ipiv, B=PhiR, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
   
   !....................................................
   ! Set MBB, MBM, and KBB from Eq. 4:
   ! Set CBB, CBM and CMM
   !....................................................
   CALL AllocAry( PhiR_T_MLL,  p%DOFR, p%DOFL, 'PhiR_T_MLL', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry( PhiR_T_CLL,  p%DOFR, p%DOFL, 'PhiR_T_MLL', ErrStat2, ErrMsg2); if(Failed()) return
      
   PhiR_T_MLL = TRANSPOSE(PhiR)
   PhiR_T_MLL = MATMUL(PhiR_T_MLL, MLL)
   MBB = MATMUL(MRL, PhiR)
   MBB = MRR + MBB + TRANSPOSE( MBB ) + MATMUL( PhiR_T_MLL, PhiR )

   PhiR_T_CLL = TRANSPOSE(PhiR)
   PhiR_T_CLL = MATMUL(PhiR_T_CLL, CLL)
   CBB = MATMUL(CRL, PhiR)
   CBB = CRR + CBB + TRANSPOSE( CBB ) + MATMUL( PhiR_T_CLL, PhiR )
   
   !MIB= MBB(p%DOFR-p%DOFI+1:p%DOFR,p%DOFR:p%DOFR) !AÑADO
      
   IF ( DOFM .EQ. 0) THEN
      MBM = 0.0_ReKi
      CBM = 0.0_ReKi
      CMM = 0.0_ReKi
   ELSE
      MBM = MATMUL( PhiR_T_MLL, PhiL(:,1:DOFM))  ! last half of operation
      MBM = MATMUL( MRL, PhiL(:,1:DOFM) ) + MBM    !This had PhiM
      
    !  MLB = MBM(p%DOFR:p%DOFR,:) !AÑADO

      CBM = MATMUL( PhiR_T_CLL, PhiL(:,1:DOFM))  ! last half of operation
      CBM = MATMUL( CRL, PhiL(:,1:DOFM) ) + CBM    !This had PhiM  

      CMM = MATMUL( MATMUL( TRANSPOSE( PhiL(:,1:DOFM) ) , CLL ) , PhiL(:,1:DOFM) )      
   ENDIF
   DEALLOCATE( PhiR_T_MLL )
   DEALLOCATE( PhiR_T_CLL )
   
   KBB = MATMUL(KRL, PhiR)   
   KBB = KBB + KRR

   IF (p%SeismicInp) THEN

   !................................
   ! set vectors related to seismic input (assumed horizontal only)
   ! The parts of the vectors that do not depend on t are set here
   !................................

   ! Set Influece vectors RR and RL

   CALL AllocAry( RL,  p%DOFL, 'Influence vector RL', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry( RR,  p%DOFR, 'Influence vector RR', ErrStat2, ErrMsg2); if(Failed()) return 
    CALL AllocAry( RI,  p%DOFI, 'Influence vector RI', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry( RB,  p%DOFR - p%DOFI, 'Influence vector RB', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry( p%RRbase,  p%DOFR - p%DOFI, 'Influence vector RRbase', ErrStat2, ErrMsg2); if(Failed()) return

   RR = 0.0_ReKi
   RI = 0.0_ReKi
   RL = 0.0_ReKi
   RB = 0.0_ReKi
   J = 0

   DO I = 1,p%DOFR
      J = J + 1
      SELECT CASE (J)
         CASE (1); RR(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)  ! The argument of COS and SIN is radians. UgDir is in º
         CASE (2); RR(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
         CASE (6); J = 0
      END SELECT 
   END DO

   p%RRbase = RR(1:p%DOFR-p%DOFI) !! RR es DOFR, p%RRbase 
   
   DO I = 1,p%DOFI
      J = J + 1
      SELECT CASE (J)
         CASE (1); RI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)  ! The argument of COS and SIN is radians. UgDir is in º
         CASE (2); RI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
         CASE (6); J = 0
      END SELECT 
   END DO

   J = 0
   DO I = 1,p%DOFL
      J = J + 1
      SELECT CASE (J)
         CASE (1); RL(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)  ! The argument of COS and SIN is radians. UgDir is in º
         CASE (2); RL(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
         CASE (6); J = 0
      END SELECT 
   END DO
   
         J = 0
       DO I = 1,p%DOFR-P%DOFI
          J = J + 1
          SELECT CASE (J)
             CASE (1); RB(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)  ! The argument of COS and SIN is radians. UgDir is in º
             CASE (2); RB(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
             CASE (6); J = 0
          END SELECT 
       END DO
       
     ! -------------------------------------------- 
     ! Set Influence vectors RR_PHI, RL_PHI, RB_PHI 
     ! --------------------------------------------
       
   CALL AllocAry( RL_PHI,  p%DOFL, 'Influence vector RL_PHI', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry( RR_PHI,  p%DOFR, 'Influence vector RR_PHI', ErrStat2, ErrMsg2); if(Failed()) return  
   CALL AllocAry( RI_PHI,  p%DOFI, 'Influence vector RI_PHI', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry( RB_PHI,  p%DOFR - p%DOFI, 'Influence vector RB_PHI', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry( p%RRbase_PHI,  p%DOFR - p%DOFI, 'Influence vector RRbase_PHI', ErrStat2, ErrMsg2); if(Failed()) return

   RR_PHI = 0.0_ReKi
   RI_PHI = 0.0_ReKi
   RL_PHI = 0.0_ReKi
   RB_PHI = 0.0_ReKi
   
   DO I = 1,p%DOFR
      J = J + 1
      SELECT CASE (J)
         CASE (1); RR_PHI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi) * (Init%Nodes(INT(p%IDR(I)/6)+1,4) - Init%Nodes(INT(p%IDR(1)/6)+1,4))  ! The argument of COS and SIN is radians. UgDir is in º
         CASE (2); RR_PHI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi) * (Init%Nodes(INT(p%IDR(I)/6)+1,4) - Init%Nodes(INT(p%IDR(1)/6)+1,4))
         CASE (4); RR_PHI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
         CASE (5); RR_PHI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)
         CASE (6); J = 0
      END SELECT 
   END DO
   
   p%RRbase_PHI = RR_PHI(1:p%DOFR-p%DOFI) !! RR_PHI es DOFR, p%RRbase_PHI es DOFI
   
   DO I = 1,p%DOFI
      J = J + 1
      SELECT CASE (J)
         CASE (1); RI_PHI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi) * (Init%Nodes(INT(p%IDI(I)/6)+1,4) - Init%Nodes(INT(p%IDR(1)/6)+1,4)) ! The argument of COS and SIN is radians. UgDir is in º
         CASE (2); RI_PHI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi) * (Init%Nodes(INT(p%IDI(I)/6)+1,4) - Init%Nodes(INT(p%IDR(1)/6)+1,4))
         CASE (4); RI_PHI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
         CASE (5); RI_PHI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)
         CASE (6); J = 0
      END SELECT 
   END DO   
      
    J = 0
      DO I = 1,p%DOFL
      J = J + 1
      SELECT CASE (J)
         CASE (1); RL_PHI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi) * (Init%Nodes(INT(p%IDL(I)/6)+1,4) - Init%Nodes(INT(p%IDR(1)/6)+1,4))! The argument of COS and SIN is radians. UgDir is in º
         CASE (2); RL_PHI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi) * (Init%Nodes(INT(p%IDL(I)/6)+1,4) - Init%Nodes(INT(p%IDR(1)/6)+1,4))
         CASE (4); RL_PHI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
         CASE (5); RL_PHI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)
         CASE (6); J = 0
      END SELECT 
   END DO
   
   J = 0
       DO I = 1,p%DOFR-p%DOFI
          J = J + 1
          SELECT CASE (J)
          CASE (1); RB_PHI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi) * (Init%Nodes(INT(p%IDR(I)/6)+1,4) - Init%Nodes(INT(p%IDR(1)/6)+1,4)) ! The argument of COS and SIN is radians. UgDir is in º
          CASE (2); RB_PHI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi) * (Init%Nodes(INT(p%IDR(I)/6)+1,4) - Init%Nodes(INT(p%IDR(1)/6)+1,4))
          CASE (4); RB_PHI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
          CASE (5); RB_PHI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)
          CASE (6); J = 0
      END SELECT 
   END DO  
             
     ! -------------------------------------------- 
     ! Set Influence vectors RR_V, RL_V, RB_V 
     ! --------------------------------------------
   
   CALL AllocAry( RL_V,  p%DOFL, 'Influence vector RL_V', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry( RR_V,  p%DOFR, 'Influence vector RR_V', ErrStat2, ErrMsg2); if(Failed()) return  
   CALL AllocAry( RI_V,  p%DOFI, 'Influence vector RI_V', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry( RB_V,  p%DOFR - p%DOFI, 'Influence vector RB_V', ErrStat2, ErrMsg2); if(Failed()) return


   RR_V = 0.0_ReKi
   RI_V = 0.0_ReKi
   RL_V = 0.0_ReKi
   RB_V = 0.0_ReKi  
   
DO I = 1,p%DOFR
      J = J + 1
      SELECT CASE (J)
         CASE (3); RR_V(I) = 1
         CASE (6); J = 0
      END SELECT 
   END DO
     
   DO I = 1,p%DOFI
      J = J + 1
      SELECT CASE (J)
         CASE (3); RI_V(I) = 1
         CASE (6); J = 0
      END SELECT 
   END DO   
      
    J = 0
      DO I = 1,p%DOFL
      J = J + 1
      SELECT CASE (J)
         CASE (3); RL_V(I) = 1
         CASE (6); J = 0
      END SELECT 
   END DO
   
   J = 0
       DO I = 1,p%DOFR-p%DOFI
          J = J + 1
          SELECT CASE (J)
          CASE (3); RB_V(I) = 1
          CASE (6); J = 0
      END SELECT 
   END DO  
        
          ! -------------------
          ! Parámetros de Sismo
          ! -------------------          
          ! Cambio RL por RM
          
      CALL AllocAry(InvPhiM, DOFM , DOFM , 'InvPhiM'      , ErrStat2, ErrMsg2); if(Failed()) return
      InvPhiM = 0.0_ReKi 
      
      IF ( p%Nmodes > 0) THEN ! Tener en cuenta invphiM=0 cuando Nmodes=0
      InvPhiM = Inv(PhiL(:,1:DOFM),DOFM) 
       
      FRIMCP = (MATMUL(CBB(p%DOFR-p%DOFI+1:p%DOFR, p%DOFR-p%DOFI+1:p%DOFR),RI) + MATMUL(CBM(p%DOFR-p%DOFI+1:p%DOFR, : ) ,matmul(InvPhiM, RL - matmul(PhiR,RR) ) ))
      
      Else 
      
      FRIMCP = (MATMUL(CBB(p%DOFR-p%DOFI+1:p%DOFR, p%DOFR-p%DOFI+1:p%DOFR),RI) + MATMUL(CBM(p%DOFR-p%DOFI+1:p%DOFR, : ) ,RL ))
         ! MATMUL( CRR + MATMUL ( TRANSPOSE(PhiR) , TRANSPOSE(CRL) ) , RR )   &
     !     + MATMUL( CRL + MATMUL ( TRANSPOSE(PhiR) , CLL ) ,  matmul(InvPhiM, RL - matmul(PhiR,RR) ) )
      
      END IF
      
                                              
      FRIMKP = MATMUL(KBB(p%DOFR-p%DOFI+1:p%DOFR, p%DOFR-p%DOFI+1:p%DOFR),RI) ! KBBb * RI
       !MATMUL( KRR + MATMUL ( TRANSPOSE(PhiR) , TRANSPOSE(KRL) ) , RR )   &
       !      + MATMUL( KRL + MATMUL ( TRANSPOSE(PhiR) , KLLcopy ) , matmul(InvPhiM, RL - matmul(PhiR,RR) ) )
          
          
      FRIMMP = MATMUL(MBB(p%DOFR-p%DOFI+1:p%DOFR,1:p%DOFR-p%DOFI),RB)  !MIB * RB
      !MATMUL(MATMUL( PhiL(:,1:DOFM), MRL+ MATMUL(TRANSPOSE(PhiR),MLL)),RB )         
   
 ! ------------------------------------------------------------------------------------------------------------------ 
 ! ESTABLEZCO FUERZAS SISMOS PARA GIRO
 
     IF ( p%Nmodes > 0) THEN ! Tener en cuenta invphiM=0 cuando Nmodes=0
       InvPhiM = Inv(PhiL(:,1:DOFM),DOFM)  
       FRIMCP_PHI = (MATMUL(CBB(p%DOFR-p%DOFI+1:p%DOFR, p%DOFR-p%DOFI+1:p%DOFR),RI_PHI) + MATMUL(CBM(p%DOFR-p%DOFI+1:p%DOFR, : ) ,matmul(InvPhiM, RL_PHI - matmul(PhiR,RR_PHI) ) )) 
     Else  
       FRIMCP_PHI = (MATMUL(CBB(p%DOFR-p%DOFI+1:p%DOFR, p%DOFR-p%DOFI+1:p%DOFR),RI_PHI) + MATMUL(CBM(p%DOFR-p%DOFI+1:p%DOFR, : ) ,RL_PHI ))
     END IF
    
   FRIMKP_PHI = MATMUL(KBB(p%DOFR-p%DOFI+1:p%DOFR, p%DOFR-p%DOFI+1:p%DOFR),RI_PHI) ! KBBb * RI
   FRIMMP_PHI = MATMUL(MBB(p%DOFR-p%DOFI+1:p%DOFR,1:p%DOFR-p%DOFI),RB_PHI)  !MIB * RB
   
   ! ------------------------------------------------------------------------------------------------------------------
   ! FUERZAS VERTICAL SISMO 
   
     IF ( p%Nmodes > 0) THEN ! Tener en cuenta invphiM=0 cuando Nmodes=0
       InvPhiM = Inv(PhiL(:,1:DOFM),DOFM)  
       FRIMCP_V = (MATMUL(CBB(p%DOFR-p%DOFI+1:p%DOFR, p%DOFR-p%DOFI+1:p%DOFR),RI_V) + MATMUL(CBM(p%DOFR-p%DOFI+1:p%DOFR, : ) ,matmul(InvPhiM, RL_V - matmul(PhiR,RR_V) ) )) 
     Else  
       FRIMCP_V = (MATMUL(CBB(p%DOFR-p%DOFI+1:p%DOFR, p%DOFR-p%DOFI+1:p%DOFR),RI_V) + MATMUL(CBM(p%DOFR-p%DOFI+1:p%DOFR, : ) ,RL_V ))
     END IF
    
   FRIMKP_V = MATMUL(KBB(p%DOFR-p%DOFI+1:p%DOFR, p%DOFR-p%DOFI+1:p%DOFR),RI_V) ! KBBb * RI
   FRIMMP_V = MATMUL(MBB(p%DOFR-p%DOFI+1:p%DOFR,1:p%DOFR-p%DOFI),RB_V)  !MIB * RB

   IF ( DOFM .EQ. 0) THEN

      FMIMKP = 0.0_ReKi
      FMIMCP = 0.0_ReKi
      FMIMMP = 0.0_ReKi
      
      FMIMKP_PHI = 0.0_ReKi
      FMIMCP_PHI = 0.0_ReKi
      FMIMMP_PHI = 0.0_ReKi 
      
      !VERTICAL
      FMIMKP_V = 0.0_ReKi
      FMIMCP_V = 0.0_ReKi
      FMIMMP_V = 0.0_ReKi 
     

   ELSE

      CALL AllocAry( OmegaDamp2 , DOFM , DOFM , 'OmegaDamp2' , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry( NOmegaM2 , DOFM , 'NomegaM2' , ErrStat2, ErrMsg2); if(Failed()) return
!      CALL AllocAry(InvPhiM, DOFM , DOFM , 'InvPhiM'      , ErrStat2, ErrMsg2); if(Failed()) return

      OmegaDamp2 = 0.0_ReKi
      NOmegaM2 = 0.0_ReKi 

      DO I = 1, DOFM
        OmegaDamp2(I,I) = 2.0_ReKi * OmegaL(I) * Init%JDampings(I) 
        NOmegaM2(I)  = -1.0_ReKi * OmegaL(I) * OmegaL(I)
      ENDDO

!      InvPhiM = Inv(PhiL(:,1:DOFM),DOFM)    !!NO SE PUEDE INVERTIR UNA MATRIZ NO CUADRADA!!
                                             !Falta ver como incluir amortiguamiento estructural   

      FMIMKP = NOmegaM2(I-1) * matmul(InvPhiM, RL - matmul(PhiR,RR) ) !KMB ES 0, SOLO APARECE KLL=OMEGA2 
      !MATMUL( TRANSPOSE( PhiL(:,1:DOFM) ) , MATMUL( TRANSPOSE(KRL) , RR ) )     &
       !         + MATMUL( TRANSPOSE( PhiL(:,1:DOFM) ) , MATMUL( KLLcopy , matmul(InvPhiM, RL - matmul(PhiR,RR) )) )

      FMIMCP = MATMUL(TRANSPOSE(CBM(p%DOFR-p%DOFI+1:p%DOFR, : )),RI) + MATMUL(CMM(:, :)+ OmegaDamp2(I-1,I-1) , MATMUL(InvPhiM, RL - MATMUL(PhiR,RR) )  )
      !MATMUL( TRANSPOSE( PhiL(:,1:DOFM) ), MATMUL( TRANSPOSE(CRL) , RR ) )      &
       !         + MATMUL( TRANSPOSE( PhiL(:,1:DOFM) ), MATMUL( CLL , matmul(InvPhiM, RL - matmul(PhiR,RR) ) ) )!+            &
       !         MATMUL( OmegaDamp2 , MATMUL( InvPhiM , RL ) )

      FMIMMP = MATMUL(TRANSPOSE(MBM(1:p%DOFR-p%DOFI,:)),RB) ! Mmb --> SAle de MMB, Mmb * RB
      !MATMUL(MATMUL(TRANSPOSE(PhiL(:,1:DOFM)),TRANSPOSE(MRL)+ MATMUL(MLL,PhiR)),RB)  !ARREGLAR FM
      
    ! -------------------------
    ! PARÁMETROS GIRO SISMO, M
    !--------------------------
    
    FMIMKP_PHI = NOmegaM2(I-1) * matmul(InvPhiM, RL_PHI - matmul(PhiR,RR_PHI) ) !KMB ES 0, SOLO APARECE KLL=OMEGA2
    
    FMIMCP_PHI = MATMUL(TRANSPOSE(CBM(p%DOFR-p%DOFI+1:p%DOFR, : )),RI_PHI) + MATMUL(CMM(:, :)+ OmegaDamp2(I-1,I-1) , MATMUL(InvPhiM, RL_PHI - MATMUL(PhiR,RR_PHI) )  )
    
    FMIMMP_PHI = MATMUL(TRANSPOSE(MBM(1:p%DOFR-p%DOFI,:)),RB_PHI) ! Mmb --> SAle de MMB, Mmb * RB
    
    ! -------------------------
    ! PARÁMETROS VERTICAL SISMO, M
    !--------------------------
    
    FMIMKP_V = NOmegaM2(I-1) * matmul(InvPhiM, RL_V - matmul(PhiR,RR_V) ) !KMB ES 0, SOLO APARECE KLL=OMEGA2
    
    FMIMCP_V = MATMUL(TRANSPOSE(CBM(p%DOFR-p%DOFI+1:p%DOFR, : )),RI_V) + MATMUL(CMM(:, :)+ OmegaDamp2(I-1,I-1) , MATMUL(InvPhiM, RL_V - MATMUL(PhiR,RR_V) )  )
    
    FMIMMP_V = MATMUL(TRANSPOSE(MBM(1:p%DOFR-p%DOFI,:)),RB_V) ! Mmb --> SAle de MMB, Mmb * RB

      DEALLOCATE(OmegaDamp2)
      DEALLOCATE(NOmegaM2)
      DEALLOCATE(InvPhiM)

   ENDIF
   
   IF (ALLOCATED(RR)) DEALLOCATE(RR)
   IF (ALLOCATED(RL)) DEALLOCATE(RL)
   IF (ALLOCATED(RB)) DEALLOCATE(RB)
   IF (ALLOCATED(RI)) DEALLOCATE(RI)
   
   IF (ALLOCATED(RR_PHI)) DEALLOCATE(RR_PHI)
   IF (ALLOCATED(RL_PHI)) DEALLOCATE(RL_PHI)
   IF (ALLOCATED(RB_PHI)) DEALLOCATE(RB_PHI)
   IF (ALLOCATED(RI_PHI)) DEALLOCATE(RI_PHI)
   
   IF (ALLOCATED(RR_V)) DEALLOCATE(RR_V)
   IF (ALLOCATED(RL_V)) DEALLOCATE(RL_V)
   IF (ALLOCATED(RB_V)) DEALLOCATE(RB_V)
   IF (ALLOCATED(RI_V)) DEALLOCATE(RI_V) 
   
   

   ENDIF  ! End IF (p%SeismicInp)

   DEALLOCATE(KLLcopy) 

CONTAINS

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CBMatrix') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   
   subroutine CleanUp()
      if (allocated(Mu        )) DEALLOCATE(Mu        )
      if (allocated(Temp      )) DEALLOCATE(Temp      )
      if (allocated(PhiR_T_MLL)) DEALLOCATE(PhiR_T_MLL)
      if (allocated(PhiR_T_CLL)) DEALLOCATE(PhiR_T_CLL)
   end subroutine
END SUBROUTINE CBMatrix

!------------------------------------------------------------------------------------------------------
!>
SUBROUTINE TrnsfTI(Init, TI, DOFI, IDI, TI2, DOFR, IDR, ErrStat, ErrMsg)
   TYPE(SD_InitType),      INTENT(IN   )  :: Init         ! Input data for initialization routine
   INTEGER(IntKi),         INTENT(IN   )  :: DOFI         ! # of DOFS of interface nodes
   INTEGER(IntKi),         INTENT(IN   )  :: DOFR         ! # of DOFS of restrained nodes (restraints and interface)
   INTEGER(IntKi),         INTENT(IN   )  :: IDI(DOFI)
   INTEGER(IntKi),         INTENT(IN   )  :: IDR(DOFR)
   REAL(ReKi),             INTENT(INOUT)  :: TI( DOFI,6)  ! matrix TI that relates the reduced matrix to the TP, 
   REAL(ReKi),             INTENT(INOUT)  :: TI2(DOFR,6)  ! matrix TI2 that relates to (0,0,0) the overall substructure mass
   INTEGER(IntKi),         INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER                                :: I, di 
   INTEGER                                :: rmndr, n
   REAL(ReKi)                             :: dx, dy, dz
   
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   TI(:,:) = 0. !Initialize     
   DO I = 1, DOFI
      di = IDI(I)
      rmndr = MOD(di, 6)
      n = CEILING(di/6.0)
      
      dx = Init%Nodes(n, 2) - Init%TP_RefPoint(1)
      dy = Init%Nodes(n, 3) - Init%TP_RefPoint(2)
      dz = Init%Nodes(n, 4) - Init%TP_RefPoint(3)
      
      SELECT CASE (rmndr)
         CASE (1); TI(I, 1:6) = (/1.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi,       dz,      -dy/)
         CASE (2); TI(I, 1:6) = (/0.0_ReKi, 1.0_ReKi, 0.0_ReKi,      -dz, 0.0_ReKi,       dx/)
         CASE (3); TI(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi,       dy,      -dx, 0.0_ReKi/)
         CASE (4); TI(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi, 0.0_ReKi/)
         CASE (5); TI(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi/)
         CASE (0); TI(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi/)
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMsg  = 'Error calculating transformation matrix TI '
            RETURN
         END SELECT
      
   ENDDO
   
   !Augment with TI2
   TI2(:,:) = 0. !Initialize 
   DO I = 1, DOFR
      di = IDR(I)
      rmndr = MOD(di, 6)
      n = CEILING(di/6.0)
      
      dx = Init%Nodes(n, 2)
      dy = Init%Nodes(n, 3) 
      dz = Init%Nodes(n, 4) 
     SELECT CASE (rmndr)
         CASE (1); TI2(I, 1:6) = (/1.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi,       dz,      -dy/)
         CASE (2); TI2(I, 1:6) = (/0.0_ReKi, 1.0_ReKi, 0.0_ReKi,      -dz, 0.0_ReKi,       dx/)
         CASE (3); TI2(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi,       dy,      -dx, 0.0_ReKi/)
         CASE (4); TI2(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi, 0.0_ReKi/)
         CASE (5); TI2(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi/)
         CASE (0); TI2(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi/)
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMsg  = 'Error calculating transformation matrix TI2 '
            RETURN
         END SELECT 
   ENDDO
   
END SUBROUTINE TrnsfTI

!------------------------------------------------------------------------------------------------------
!> Return eigenvalues, Omega, and eigenvectors, Phi, 
SUBROUTINE EigenSolve(K, M, nDOF, NOmega, Reduced, Init,p, Phi, Omega, ErrStat, ErrMsg )
   USE NWTC_ScaLAPACK, only: ScaLAPACK_LASRT
   INTEGER,                INTENT(IN   )    :: nDOF                               ! Total degrees of freedom of the incoming system
   REAL(ReKi),             INTENT(IN   )    :: K(nDOF, nDOF)                      ! stiffness matrix 
   REAL(ReKi),             INTENT(IN   )    :: M(nDOF, nDOF)                      ! mass matrix 
   INTEGER,                INTENT(IN   )    :: NOmega                             ! RRD: no. of requested eigenvalues
   LOGICAL,                INTENT(IN   )    :: Reduced                            ! Whether or not to reduce matrices, this will be removed altogether later, when reduction will be done apriori
   TYPE(SD_InitType),      INTENT(IN   )    :: Init  
   TYPE(SD_ParameterType), INTENT(IN   )    :: p  
   REAL(ReKi),             INTENT(  OUT)    :: Phi(nDOF, NOmega)                  ! RRD: Returned Eigenvectors
   REAL(ReKi),             INTENT(  OUT)    :: Omega(NOmega)                      ! RRD: Returned Eigenvalues
   INTEGER(IntKi),         INTENT(  OUT)    :: ErrStat                            ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT)    :: ErrMsg                             ! Error message if ErrStat /= ErrID_None
   ! LOCALS         
   REAL(LAKi), ALLOCATABLE                   :: Omega2(:)                         !RRD: Eigen-values new system
! note: SGGEV seems to have memory issues in certain cases. The eigenvalues seem to be okay, but the eigenvectors vary wildly with different compiling options.
!       DGGEV seems to work better, so I'm making these variables LAKi (which is set to R8Ki for now)   - bjj 4/25/2014
   REAL(LAKi), ALLOCATABLE                   :: Kred(:,:), Mred(:,:) 
   REAL(LAKi), ALLOCATABLE                   :: WORK (:),  VL(:,:), VR(:,:), ALPHAR(:), ALPHAI(:), BETA(:) ! eigensolver variables
   INTEGER                                   :: i  
   INTEGER                                   :: N, LWORK                          !variables for the eigensolver
   INTEGER,    ALLOCATABLE                   :: KEY(:)
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
      
   ErrStat = ErrID_None
   ErrMsg  = ''
         
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
   IF (Reduced) THEN !bjj: i.e., We need to reduce; it's not reduced yet
      ! First I need to remove constrained nodes DOFs
      ! This is actually done when we are printing out the 'full' set of eigenvalues
      CALL ReduceKMdofs(Kred,K,nDOF, Init,p, ErrStat2, ErrMsg2 ); if(Failed()) return
      CALL ReduceKMdofs(Mred,M,nDOF, Init,p, ErrStat2, ErrMsg2 ); if(Failed()) return
      N=SIZE(Kred,1)   
   ELSE
      ! This is actually done whe we are generating the CB-reduced set of eigenvalues, so the the variable 'Reduced' can be a bit confusing. GJH 8/1/13
      N=SIZE(K,1)
      CALL AllocAry( Kred, n, n, 'Kred', ErrStat2, ErrMsg2 ); if(Failed()) return
      CALL AllocAry( Mred, n, n, 'Mred', ErrStat2, ErrMsg2 ); if(Failed()) return
      Kred=REAL( K, LAKi )
      Mred=REAL( M, LAKi )
   ENDIF
   ! Note:  NOmega must be <= N, which is the length of Omega2, Phi!
   IF ( NOmega > N ) THEN
      CALL SetErrStat(ErrID_Fatal,"NOmega must be less than or equal to N",ErrStat,ErrMsg,'EigenSolve')
      CALL CleanupEigen()
      RETURN
   END IF

   ! allocate working arrays and return arrays for the eigensolver
   LWORK=8*N + 16  !this is what the eigensolver wants  >> bjj: +16 because of MKL ?ggev documenation ( "lwork >= max(1, 8n+16) for real flavors"), though LAPACK documenation says 8n is fine
   !bjj: there seems to be a memory problem in *GGEV, so I'm making the WORK array larger to see if I can figure it out
   CALL AllocAry( Work,   lwork,     'Work',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve') 
   CALL AllocAry( Omega2, n,         'Omega2', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( ALPHAR, n,         'ALPHAR', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( ALPHAI, n,         'ALPHAI', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( BETA,   n,         'BETA',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( VR,     n,  n,     'VR',     ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( VL,     n,  n,     'VR',     ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'EigenSolve')
   CALL AllocAry( KEY,    n,         'KEY',    ErrStat2, ErrMsg2 ); if(Failed()) return
    
   CALL  LAPACK_ggev('N','V',N ,Kred, Mred, ALPHAR, ALPHAI, BETA, VL, VR, work, lwork, ErrStat2, ErrMsg2)
   if(Failed()) return
   !if (.not. reduced) call wrmatrix(REAL(VR,ReKi),77,'ES15.8e2')    
   ! bjj: This comes from the LAPACK documentation:
   !   Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) may easily over- or underflow, and BETA(j) may even be zero.
   !   Thus, the user should avoid naively computing the ratio alpha/beta.  However, ALPHAR and ALPHAI will be always less
   !   than and usually comparable with norm(A) in magnitude, and BETA always less than and usually comparable with norm(B).    
   ! Omega2=ALPHAR/BETA  !Note this may not be correct if ALPHAI<>0 and/or BETA=0 TO INCLUDE ERROR CHECK, also they need to be sorted
   DO I=1,N !Initialize the key and calculate Omega2
      KEY(I)=I
      IF ( EqualRealNos(Beta(I),0.0_LAKi) ) THEN
         Omega2(I) = HUGE(Omega2)  ! bjj: should this be an error?
      ELSE
         Omega2(I) = REAL( ALPHAR(I)/BETA(I), ReKi )
      END IF           
   ENDDO  
   CALL ScaLAPACK_LASRT('I',N,Omega2,key,ErrStat2,ErrMsg2); if(Failed()) return
    
   !we need to rearrange eigenvectors based on sorting of Omega2
   !Now rearrange VR based on the new key, also I might have to scale the eigenvectors following generalized mass =idnetity criterion, also if i reduced the matrix I will need to re-expand the eigenvector
   ! ALLOCATE(normcoeff(N,N), STAT = ErrStat )
   ! result1 = matmul(Mred2,VR)
   ! result2 = matmul(transpose(VR),result1)
   ! normcoeff=sqrt(result2)  !This should be a diagonal matrix which contains the normalization factors
   !normcoeff=sqrt(matmul(transpose(VR),matmul(Mred2,VR)))  !This should be a diagonal matrix which contains the normalization factors
   VL=VR  !temporary storage for sorting VR
   DO I=1,N 
      !VR(:,I)=VL(:,KEY(I))/normcoeff(KEY(I),KEY(I))  !reordered and normalized
      VR(:,I)=VL(:,KEY(I))  !just reordered as Huimin had a normalization outside of this one
   ENDDO
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

   ! --- Finish EigenSolve
   ! Note:  NOmega must be <= N, which is the length of Omega2, Phi!
   Omega=SQRT( Omega2(1:NOmega) ) !Assign my new Omega and below my new Phi (eigenvectors) [eigenvalues are actually the square of omega]
   IF ( Reduced ) THEN ! this is called for the full system Eigenvalues:
      !Need to expand eigenvectors for removed DOFs, setting Phi 
      CALL UnReduceVRdofs(VR(:,1:NOmega),Phi,N,NOmega, Init,p, ErrStat2, ErrMsg2 ) ; if(Failed()) return
   ELSE ! IF (.NOT.(Reduced)) THEN !For the time being Phi gets updated only when CB eigensolver is requested. I need to fix it for the other case (full fem) and then get rid of the other eigensolver, this implies "unreducing" the VR
       ! This is done as part of the CB-reduced eigensolve
      Phi=REAL( VR(:,1:NOmega), ReKi )   ! eigenvectors
   ENDIF  
   
   CALL CleanupEigen()
   RETURN

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'EigenSolve') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUpEigen()
   END FUNCTION Failed

   SUBROUTINE CleanupEigen()
      IF (ALLOCATED(Work)  ) DEALLOCATE(Work)
      IF (ALLOCATED(Omega2)) DEALLOCATE(Omega2)  !bjj: break in Debug_Doub
      IF (ALLOCATED(ALPHAR)) DEALLOCATE(ALPHAR)
      IF (ALLOCATED(ALPHAI)) DEALLOCATE(ALPHAI)
      IF (ALLOCATED(BETA)  ) DEALLOCATE(BETA)
      IF (ALLOCATED(VR)    ) DEALLOCATE(VR)
      IF (ALLOCATED(VL)    ) DEALLOCATE(VL)
      IF (ALLOCATED(KEY)   ) DEALLOCATE(KEY)
      IF (ALLOCATED(Kred)  ) DEALLOCATE(Kred)
      IF (ALLOCATED(Mred)  ) DEALLOCATE(Mred)
   END SUBROUTINE CleanupEigen
  
END SUBROUTINE EigenSolve

!------------------------------------------------------------------------------------------------------
!> Calculate Kred from K after removing consstrained node DOFs from the full M and K matrices
!!Note it works for constrained nodes, still to see how to make it work for interface nodes if needed
SUBROUTINE ReduceKMdofs(Kred,K,TDOF, Init,p, ErrStat, ErrMsg )
   TYPE(SD_InitType),      INTENT(  in)  :: Init  
   TYPE(SD_ParameterType), INTENT(  in)  :: p  
   INTEGER,                INTENT(IN   ) :: TDOF           ! Size of matrix K (total DOFs)                              
   REAL(ReKi),             INTENT(IN   ) :: K(TDOF, TDOF)  ! full matrix
   REAL(LAKi),ALLOCATABLE, INTENT(  OUT) :: Kred(:,:)      ! reduced matrix
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat        ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   !locals
   INTEGER                               :: I, J           ! counters into full or reduced matrix
   INTEGER                               :: L              ! number of DOFs to eliminate 
   INTEGER, ALLOCATABLE                  :: idx(:)         ! vector to map reduced matrix to full matrix
   INTEGER                               :: NReactDOFs
   INTEGER                               :: DOF_reduced
   INTEGER                               :: ErrStat2
   CHARACTER(1024)                       :: ErrMsg2
   
   ErrStat = ErrID_None
   ErrMsg  = ''    
  
   NReactDOFs = p%NReact*6 !p%DOFC
   IF (NReactDOFs > TDOF) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = 'ReduceKMdofs:invalid matrix sizes.'
      RETURN
   END IF
  
   CALL AllocAry(idx,  TDOF, 'idx',  ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,'ReduceKMdofs')
   IF (ErrStat >= AbortErrLev) THEN
      RETURN
   END IF   
   
   ! Calculate how many rows/columns need to be eliminated:
   DO I = 1, TDOF       
      idx(I) = I
   END DO
         
   L = 0
   DO I = 1, NReactDOFs  !Cycle on reaction DOFs      
      IF (Init%BCs(I, 2) == 1) THEN
         L=L+1 !number of DOFs to eliminate
         idx( Init%BCs(I, 1) ) = 0 ! Eliminate this one
      END IF    
   END DO   
   
   ! Allocate the output matrix and the index mapping array
   DOF_reduced = TDOF-L   
   CALL AllocAry(Kred, DOF_reduced, DOF_reduced, 'Kred', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,'ReduceKMdofs')
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp()
      RETURN
   END IF
   
   ! set the indices we want to keep (i.e., a mapping from reduced to full matrix)
   J  = 1   
   DO I=1,TDOF
      idx(J) = idx(I)      
      IF ( idx(J) /= 0 ) J = J + 1  
       
   END DO
         
   ! Remove rows and columns from every row/column in full matrix where Init%BC(:,2) == 1,
   ! using the mapping created above. (This is a symmetric matrix.)
   DO J = 1, DOF_reduced  !Cycle on reaction DOFs      
      DO I = 1, DOF_reduced  !Cycle on reaction DOFs      
         Kred(I,J) = REAL( K( idx(I), idx(J) ), LAKi )
      END DO
   END DO 

   ! clean up local variables:
   CALL CleanUp()
CONTAINS
   subroutine CleanUp()
      IF (ALLOCATED(idx)) DEALLOCATE(idx)
   end subroutine
END SUBROUTINE ReduceKMdofs

!------------------------------------------------------------------------------------------------------
!> Augments VRred to VR for the constrained DOFs, somehow reversing what ReducedKM did for matrices
!Note it works for constrained nodes, still to see how to make it work for interface nodes if needed
SUBROUTINE UnReduceVRdofs(VRred,VR,rDOF,rModes, Init,p, ErrStat, ErrMsg )
   TYPE(SD_InitType),      INTENT(in   ) :: Init  
   TYPE(SD_ParameterType), INTENT(in   ) :: p  
   INTEGER,                INTENT(IN   ) :: rDOF ,RModes  !retained DOFs after removing restrained DOFs and retained modes 
   REAL(LAKi),             INTENT(IN   ) :: VRred(rDOF, rModes)  !eigenvector matrix with restrained DOFs removed
   REAL(ReKi),             INTENT(INOUT) :: VR(:,:) !eigenvalues including the previously removed DOFs
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   !locals
   INTEGER,   ALLOCATABLE   :: idx(:)
   INTEGER                  :: I, I2, L  !counters; I,I2 should be long, L short

   ErrStat = ErrID_None
   ErrMsg  = ''    
  
   ALLOCATE(idx(p%NReact*6), STAT = ErrStat )  !it contains row/col index that was originally eliminated when applying restraints
   idx=0 !initialize
   L=0 !initialize
   DO I = 1, p%NReact*6  !Cycle on reaction DOFs
       IF (Init%BCs(I, 2) == 1) THEN
           idx(I)=Init%BCs(I, 1) !row/col index that was originally eliminated when applying restraints
           L=L+1 !number of DOFs to eliminate
       ENDIF    
   ENDDO
!  PRINT *, '    rDOF+L=',rDOF+L, 'SIZE(Phi2)=',SIZE(VR,1)
!  ALLOCATE(VR(rDOF+L,rModes), STAT = ErrStat )  !Restored eigenvectors with restrained node DOFs included
   VR=0.!Initialize

   I2=1 !Initialize 
   DO I=1,rDOF+L  !This loop inserts Vred in VR in all but the removed DOFs
      IF (ALL((idx-I).NE.0)) THEN
         VR(I,:)=REAL( VRred(I2,:), ReKi ) ! potentially change of precision
         I2=I2+1  !Note this counter gets updated only if we insert Vred rows into VR
      ENDIF   
   ENDDO
END SUBROUTINE UnReduceVRdofs

!------------------------------------------------------------------------------------------------------
SUBROUTINE CBApplyConstr(DOFI, DOFR, DOFM,  DOFL,  &
                         MBB , MBM , KBB ,                      &
                         CBB , CBM , CMM , PHiR , FGR ,         &
                         FRIMKP , FRIMCP , FRIMMP,              & !AÑADO FM
                         FMIMKP , FMIMCP , FMIMMP,              & !AÑADO FM
                         FRIMKPb, FRIMCPb, FRIMMPb,             & !AÑADO FM
                         FMIMKPb, FMIMCPb, FMIMMPb,             & !AÑADO FM
                         !AÑADO GIROS
                         FRIMKP_PHI , FRIMCP_PHI , FRIMMP_PHI,  & 
                         FMIMKP_PHI , FMIMCP_PHI , FMIMMP_PHI,  & 
                         FRIMKP_PHIb, FRIMCP_PHIb, FRIMMP_PHIb, & 
                         FMIMKP_PHIb, FMIMCP_PHIb, FMIMMP_PHIb, &  
                         ! AÑADO VERTICAL
                         FRIMKP_V , FRIMCP_V , FRIMMP_V,        & 
                         FMIMKP_V , FMIMCP_V , FMIMMP_V,        & 
                         FRIMKP_Vb, FRIMCP_Vb, FRIMMP_Vb,       & 
                         FMIMKP_Vb, FMIMCP_Vb, FMIMMP_Vb,       &                         
                         MBBb, MBMb, KBBb,                      &
                         CBBb, CBMb, CMMb, PHiRb, PhiRbase, FGRb, p)

   TYPE(SD_ParameterType), INTENT(INOUT)  :: p           ! Parameters 
   INTEGER(IntKi),         INTENT(IN   )  ::  DOFR, DOFI, DOFM, DOFL
   REAL(ReKi),             INTENT(IN   )  ::  FGR(DOFR)
   REAL(ReKi),             INTENT(IN   )  ::  FRIMKP(DOFI)
   REAL(ReKi),             INTENT(IN   )  ::  FMIMKP(DOFM)
   REAL(ReKi),             INTENT(IN   )  ::  FRIMCP(DOFI)
   REAL(ReKi),             INTENT(IN   )  ::  FMIMCP(DOFM)
   REAL(ReKi),             INTENT(IN   )  ::  FRIMMP(DOFI) ! AÑADO FM
   REAL(ReKi),             INTENT(IN   )  ::  FMIMMP(DOFM) ! AÑADO FM 
   !AÑADO GIROS
   REAL(ReKi),             INTENT(IN   )  ::  FRIMKP_PHI(DOFI)
   REAL(ReKi),             INTENT(IN   )  ::  FMIMKP_PHI(DOFM)
   REAL(ReKi),             INTENT(IN   )  ::  FRIMCP_PHI(DOFI)
   REAL(ReKi),             INTENT(IN   )  ::  FMIMCP_PHI(DOFM)
   REAL(ReKi),             INTENT(IN   )  ::  FRIMMP_PHI(DOFI) 
   REAL(ReKi),             INTENT(IN   )  ::  FMIMMP_PHI(DOFM)  
   !AÑADO VERTICAL
   REAL(ReKi),             INTENT(IN   )  ::  FRIMKP_V(DOFI)
   REAL(ReKi),             INTENT(IN   )  ::  FMIMKP_V(DOFM)
   REAL(ReKi),             INTENT(IN   )  ::  FRIMCP_V(DOFI)
   REAL(ReKi),             INTENT(IN   )  ::  FMIMCP_V(DOFM)
   REAL(ReKi),             INTENT(IN   )  ::  FRIMMP_V(DOFI) 
   REAL(ReKi),             INTENT(IN   )  ::  FMIMMP_V(DOFM)     
   REAL(ReKi),             INTENT(IN   )  ::  MBB(DOFR, DOFR)
   REAL(ReKi),             INTENT(IN   )  ::  MBM(DOFR, DOFM)
   REAL(ReKi),             INTENT(IN   )  ::  KBB(DOFR, DOFR)
   REAL(ReKi),             INTENT(IN   )  ::  CBB(DOFR, DOFR)
   REAL(ReKi),             INTENT(IN   )  ::  CBM(DOFR, DOFM)
   REAL(ReKi),             INTENT(IN   )  ::  CMM(DOFM, DOFM)
   REAL(ReKi),             INTENT(IN   )  :: PhiR(DOFL, DOFR)   
   REAL(ReKi),             INTENT(  OUT)  ::  MBBb(DOFI, DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  KBBb(DOFI, DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  MBMb(DOFI, DOFM)
   REAL(ReKi),             INTENT(  OUT)  ::  CBBb(DOFI, DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  CBMb(DOFI, DOFM)
   REAL(ReKi),             INTENT(  OUT)  ::  CMMb(DOFM, DOFM)
   REAL(ReKi),             INTENT(  OUT)  ::  FGRb(DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  FRIMKPb(DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  FMIMKPb(DOFM)
   REAL(ReKi),             INTENT(  OUT)  ::  FRIMCPb(DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  FMIMCPb(DOFM)
   REAL(ReKi),             INTENT(  OUT)  ::  FRIMMPb(DOFI) !AÑADO FM
   REAL(ReKi),             INTENT(  OUT)  ::  FMIMMPb(DOFM) !AÑADO FM  
   !AÑADO GIROS
   REAL(ReKi),             INTENT(  OUT)  ::  FRIMKP_PHIb(DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  FMIMKP_PHIb(DOFM)
   REAL(ReKi),             INTENT(  OUT)  ::  FRIMCP_PHIb(DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  FMIMCP_PHIb(DOFM)
   REAL(ReKi),             INTENT(  OUT)  ::  FRIMMP_PHIb(DOFI) 
   REAL(ReKi),             INTENT(  OUT)  ::  FMIMMP_PHIb(DOFM)   
   !AÑADO VERTICALb
   REAL(ReKi),             INTENT(  OUT)  ::  FRIMKP_Vb(DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  FMIMKP_Vb(DOFM)
   REAL(ReKi),             INTENT(  OUT)  ::  FRIMCP_Vb(DOFI)
   REAL(ReKi),             INTENT(  OUT)  ::  FMIMCP_Vb(DOFM)
   REAL(ReKi),             INTENT(  OUT)  ::  FRIMMP_Vb(DOFI) 
   REAL(ReKi),             INTENT(  OUT)  ::  FMIMMP_Vb(DOFM)    
   REAL(ReKi),             INTENT(  OUT)  ::  PhiRb(DOFL, DOFI)   
   REAL(ReKi),             INTENT(  OUT)  ::  PhiRbase(DOFL, DOFR-DOFI) 
      
   MBBb  = MBB(DOFR-DOFI+1:DOFR, DOFR-DOFI+1:DOFR) 
   KBBb  = KBB(DOFR-DOFI+1:DOFR, DOFR-DOFI+1:DOFR)  
   CBBb  = CBB(DOFR-DOFI+1:DOFR, DOFR-DOFI+1:DOFR)   
IF (DOFM > 0) THEN   
   MBMb  = MBM(DOFR-DOFI+1:DOFR, :               )
   CBMb  = CBM(DOFR-DOFI+1:DOFR, :               )
   CMMb  = CMM(:, :)
   FMIMKPb = FMIMKP(:)
   FMIMCPb = FMIMCP(:)
   FMIMMPb = FMIMMP(:) !AÑADO FM
 ! AÑADO GIRO ---------------------------------------
   FMIMKP_PHIb = FMIMKP_PHI(:)
   FMIMCP_PHIb = FMIMCP_PHI(:)
   FMIMMP_PHIb = FMIMMP_PHI(:)
 ! AÑADO VERTICAL -----------------------------------
   FMIMKP_Vb = FMIMKP_V(:)
   FMIMCP_Vb = FMIMCP_V(:)
   FMIMMP_Vb = FMIMMP_V(:)   
  
END IF
   FGRb  = FGR(DOFR-DOFI+1:DOFR )
   FRIMKPb  = FRIMKP
   !FRIMKP(DOFR-DOFI+1:DOFR )
   FRIMCPb  = FRIMCP
   !FRIMCP(DOFR-DOFI+1:DOFR )
   FRIMMPb  = FRIMMP
   !FRIMMP(DOFR-DOFI+1:DOFR ) !AÑADO FM
   
   ! AÑADO GIRO --------------------------------------
   FRIMKP_PHIb  = FRIMKP_PHI
   FRIMCP_PHIb  = FRIMCP_PHI
   FRIMMP_PHIb  = FRIMMP_PHI
   ! AÑADO VERTICAL-----------------------------------
   FRIMKP_Vb  = FRIMKP_V
   FRIMCP_Vb  = FRIMCP_V
   FRIMMP_Vb  = FRIMMP_V
   
   PhiRb = PhiR(              :, DOFR-DOFI+1:DOFR)

   IF (p%SeismicInp) THEN
      PhiRbase = PhiR(        :,1:DOFR-DOFI)
   ENDIF
   
END SUBROUTINE CBApplyConstr

!------------------------------------------------------------------------------------------------------
SUBROUTINE SetParameters(Init, p, MBBb, MBmb, KBBb, CBBb, CBmb, CMMb, FGRb, FRIMKPb, FRIMCPb, FRIMMPb , FMIMKPb, FMIMCPb, FMIMMPb, &
           FRIMKP_PHIb, FRIMCP_PHIb, FRIMMP_PHIb , FMIMKP_PHIb, FMIMCP_PHIb, FMIMMP_PHIb,                                          &
           FRIMKP_Vb, FRIMCP_Vb, FRIMMP_Vb , FMIMKP_Vb, FMIMCP_Vb, FMIMMP_Vb, PhiRb, PhiRbase, OmegaL, FGL, PhiL, ErrStat, ErrMsg) !AÑADO FM 
! AÑADO GIRO, AÑADO VERTICAL 

   TYPE(SD_InitType),        INTENT(IN   )   :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType),   INTENT(INOUT)   :: p            ! Parameters
   REAL(ReKi),               INTENT(IN   )   :: MBBb(  p%DOFI, p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: MBMb(  p%DOFI, p%Nmodes)
   REAL(ReKi),               INTENT(IN   )   :: KBBb(  p%DOFI, p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: CBBb(  p%DOFI, p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: CBMb(  p%DOFI, p%Nmodes)
   REAL(ReKi),               INTENT(IN   )   :: CMMb(  p%Nmodes, p%Nmodes)
   REAL(ReKi),               INTENT(IN   )   :: PhiL ( p%DOFL, p%DOFL)   
   REAL(ReKi),               INTENT(IN   )   :: PhiRb( p%DOFL, p%DOFI) 
   REAL(ReKi),               INTENT(IN   )   :: PhiRbase( p%DOFL, p%DOFR-p%DOFI)   
   REAL(ReKi),               INTENT(IN   )   :: OmegaL(p%DOFL)   
   REAL(ReKi),               INTENT(IN   )   :: FGRb(p%DOFI) 
   REAL(ReKi),               INTENT(IN   )   :: FRIMKPb(p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: FRIMCPb(p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: FRIMMPb(p%DOFI) !AÑADO FM 
   REAL(ReKi),               INTENT(IN   )   :: FMIMKPb(p%Nmodes)
   REAL(ReKi),               INTENT(IN   )   :: FMIMCPb(p%Nmodes)
   REAL(ReKi),               INTENT(IN   )   :: FMIMMPb(p%Nmodes) !AÑADO FM  
   ! AÑADO GIRO
   REAL(ReKi),               INTENT(IN   )   :: FRIMKP_PHIb(p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: FRIMCP_PHIb(p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: FRIMMP_PHIb(p%DOFI) 
   REAL(ReKi),               INTENT(IN   )   :: FMIMKP_PHIb(p%Nmodes)
   REAL(ReKi),               INTENT(IN   )   :: FMIMCP_PHIb(p%Nmodes)
   REAL(ReKi),               INTENT(IN   )   :: FMIMMP_PHIb(p%Nmodes) 
   ! AÑADO VERTICAL
   REAL(ReKi),               INTENT(IN   )   :: FRIMKP_Vb(p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: FRIMCP_Vb(p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: FRIMMP_Vb(p%DOFI) 
   REAL(ReKi),               INTENT(IN   )   :: FMIMKP_Vb(p%Nmodes)
   REAL(ReKi),               INTENT(IN   )   :: FMIMCP_Vb(p%Nmodes)
   REAL(ReKi),               INTENT(IN   )   :: FMIMMP_Vb(p%Nmodes)           
   REAL(ReKi),               INTENT(IN   )   :: FGL(p%DOFL)
   INTEGER(IntKi),           INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   REAL(ReKi)                                :: TI_transpose(TPdofL,p%DOFI) !bjj: added this so we don't have to take the transpose 5+ times
   INTEGER(IntKi)                            :: I
   integer(IntKi)                            :: n                          ! size of jacobian in AM2 calculation
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   CHARACTER(*), PARAMETER                   :: RoutineName = 'SetParameters'
   
   ErrStat = ErrID_None 
   ErrMsg  = ''
      
   TI_transpose =  TRANSPOSE(p%TI) 

   ! Store FGL for later processes
   IF (p%SttcSolve) THEN     
       p%FGL = FGL  
   ENDIF     
      
   ! block element of D2 matrix (D2_21, D2_42, & part of D2_62)
   p%PhiRb_TI = MATMUL(PhiRb, p%TI)
   
   !...............................
   ! equation 46-47 (used to be 9):
   !...............................
   p%MBB = MATMUL( MATMUL( TI_transpose, MBBb ), p%TI) != MBBt
   p%KBB = MATMUL( MATMUL( TI_transpose, KBBb ), p%TI) != KBBt
   p%CBB = MATMUL( MATMUL( TI_transpose, CBBb ), p%TI) != CBBt

   !p%D1_15=-TI_transpose  !this is 6x6NIN
   IF ( p%NModes > 0 ) THEN ! These values don't exist for DOFM=0; i.e., p%NModes == 0
         ! p%MBM = MATMUL( TRANSPOSE(p%TI), MBmb )    != MBMt
      CALL LAPACK_gemm( 'T', 'N', 1.0_ReKi, p%TI, MBMb, 0.0_ReKi, p%MBM, ErrStat2, ErrMsg2) != MBMt
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%MBM')
      CALL LAPACK_gemm( 'T', 'N', 1.0_ReKi, p%TI, CBMb, 0.0_ReKi, p%CBM, ErrStat2, ErrMsg2) != CBMt
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%CBM')

      p%CMM = CMMb  ! No transformation needed for MM 
      
      p%MMB = TRANSPOSE( p%MBM )                          != MMBt
      p%CMB = TRANSPOSE( p%CBM )                          != CMBt

    ! --------------------------------------------------------------
    ! FORMULATION OF STATE-SPACE FORMULATION MATRICES AND VECTORS 
    ! -------------------------------------------------------------

      p%PhiM  = PhiL(:,1:p%Nmodes)
      
      ! A_21, A_22 (these are diagonal matrices. bjj: I am storing them as arrays instead of full matrices)
      p%NOmegaM2      = -1.0_ReKi * OmegaL(1:p%Nmodes) * OmegaL(1:p%Nmodes)          ! OmegaM is a one-dimensional array
      p%N2OmegaMJDamp = -2.0_ReKi * OmegaL(1:p%Nmodes) * Init%JDampings(1:p%Nmodes)  ! Init%JDampings is also a one-dimensional array
   
      ! B_23, B_24
      !p%PhiM_T =  TRANSPOSE( p%PhiM  )
   
      ! FX
      ! p%FX = MATMUL( p%PhiM_T, FGL ) != MATMUL( TRANSPOSE(PhiM), FGL )
      p%FX = MATMUL( FGL, p%PhiM ) != MATMUL( TRANSPOSE(PhiM), FGL ) because FGL is 1-D
   
      ! C1_11, C1_12  ( see eq 15 [multiply columns by diagonal matrix entries for diagonal multiply on the left])   
      DO I = 1, p%Nmodes ! if (p%NModes=p%qmL=DOFM == 0), this loop is skipped
         p%C1_11(:, I) = p%MBM(:, I)*p%NOmegaM2(I)              
         p%C1_12(:, I) = p%MBM(:, I)*p%N2OmegaMJDamp(I)  
      ENDDO   

      !IF ( p%Nmodes > 0) THEN
      CALL LAPACK_GEMM( 'N', 'N', -1.0_ReKi, p%MBM,   p%CMM,  1.0_ReKi, p%C1_12, ErrStat2, ErrMsg2 )  ! p%D1_12 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      p%C1_12 = p%C1_12 + p%CBM
      !ENDIF
   
      ! D1_12, D1_13, D1_14 (with retained modes)

      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%MBM,   p%CBM,  0.0_ReKi, p%D1_12, ErrStat2, ErrMsg2 )  ! p%D1_12 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      p%D1_12 = p%CBB - p%D1_12

      !p%D1_13 = p%MBB - MATMUL( p%MBM, p%MMB )
      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%MBM,   p%MBM,  0.0_ReKi, p%D1_13, ErrStat2, ErrMsg2 )  ! p%D1_13 = MATMUL( p%MBM, p%MMB )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      p%D1_13 = p%MBB - p%D1_13

      !p%D1_14 = MATMUL( p%MBM, p%PhiM_T ) - MATMUL( TI_transpose, TRANSPOSE(PHiRb))  
      
     
      CALL LAPACK_GEMM( 'T', 'T', 1.0_ReKi, p%TI,   PHiRb,  0.0_ReKi, p%D1_14, ErrStat2, ErrMsg2 )  ! p%D1_14 = MATMUL( TRANSPOSE(TI), TRANSPOSE(PHiRb))  
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%MBM, p%PhiM, -1.0_ReKi, p%D1_14, ErrStat2, ErrMsg2 )  ! p%D1_14 = MATMUL( p%MBM, TRANSPOSE(p%PhiM) ) - p%D1_14 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   
      ! FY (with retained modes)
      p%FY =    MATMUL( p%MBM, p%FX ) &  
              - MATMUL( TI_transpose, ( FGRb + MATMUL( TRANSPOSE(PhiRb), FGL) ) ) 
      
      ! C2_21, C2_42
      ! C2_61, C2_62
      DO I = 1, p%Nmodes ! if (p%NModes=p%qmL=DOFM == 0), this loop is skipped
         p%C2_61(:, i) = p%PhiM(:, i)*p%NOmegaM2(i)
         p%C2_62(:, i) = p%PhiM(:, i)*p%N2OmegaMJDamp(i)
      ENDDO   

      !IF ( p%Nmodes > 0) 
      p%C2_62 = p%C2_62  - MATMUL( p%PhiM , p%CMM )
      !ENDIF
      
      ! D2_63, D2_63, D2_64 
      p%D2_62 = MATMUL( p%PhiM, p%CMB ) ! ¿Debería ser negativo?
      p%D2_62 = - p%D2_62
      p%D2_63 = MATMUL( p%PhiM, p%MMB )
      p%D2_63 = p%PhiRb_TI - p%D2_63

      !p%D2_64 = MATMUL( p%PhiM, p%PhiM_T )  !bjj: why does this use stack space?
      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%PhiM, p%PhiM, 0.0_ReKi, p%D2_64, ErrStat2, ErrMsg2 ) !bjj: replaced MATMUL with this routine to avoid issues with stack size
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      ! F2_61
      p%F2_61 = MATMUL( p%D2_64, FGL )       
                              
     !Now calculate a Jacobian used when AM2 is called and store in parameters    
      IF (p%IntMethod .EQ. 4) THEN       ! Allocate Jacobian if AM2 is requested & if there are states (p%qmL > 0)
         n=2*p%qmL
         CALL AllocAry( p%AM2Jac, n, n, 'p%AM2InvJac', ErrStat2, ErrMsg2 ); if(Failed()) return
         CALL AllocAry( p%AM2JacPiv, n, 'p%AM2JacPiv', ErrStat2, ErrMsg2 ); if(Failed()) return
         
         ! First we calculate the Jacobian:
         ! (note the Jacobian is first stored as p%AM2InvJac)
         p%AM2Jac=0.
         DO i=1,p%qmL
            p%AM2Jac(i+p%qmL,i      )=p%SDdeltaT/2.*p%NOmegaM2(i)      !J21   
            p%AM2Jac(i+p%qmL,i+p%qmL)=p%SDdeltaT/2.*p%N2OmegaMJDamp(i) !J22 -initialize
         END DO
      
         DO I=1,p%qmL
            p%AM2Jac(I,I)=-1.  !J11
            p%AM2Jac(I,p%qmL+I)=p%SDdeltaT/2.  !J12
            p%AM2Jac(p%qmL+I,p%qmL+I)=p%AM2Jac(p%qmL+I,p%qmL+I)-1  !J22 complete
         ENDDO
         ! Now need to factor it:        
         !I think it could be improved and made more efficient if we can say the matrix is positive definite
         CALL LAPACK_getrf( n, n, p%AM2Jac, p%AM2JacPiv, ErrStat2, ErrMsg2); if(Failed()) return
      END IF     
      
   ELSE ! no retained modes, so 
      ! OmegaM, JDampings, PhiM, MBM, MMB, FX , x don't exist in this case
      ! p%F2_61, p%D2_64 are zero in this case so we simplify the equations in the code, omitting these variables
      ! p%D2_63 = p%PhiRb_TI in this case so we simplify the equations in the code, omitting storage of this variable
      ! p%D1_12 = p%CBB in this case so we simplify the equations in the code, omitting storage of this variable
      ! p%D1_13 = p%MBB in this case so we simplify the equations in the code, omitting storage of this variable
      
      ! D1_14 (with 0 retained modes)
      p%D1_14 = - MATMUL( TI_transpose, TRANSPOSE(PHiRb))  

      ! FY (with 0 retained modes)
      p%FY    = - MATMUL( TI_transpose, ( FGRb + MATMUL( TRANSPOSE(PhiRb), FGL) ) ) 
                  
   END IF

   IF (p%SeismicInp) THEN
    p%FRIMKP = FRIMKPb
    p%FRIMCP = FRIMCPb
    p%FRIMMP = FRIMMPb
    ! AÑADO GIRO --------------------------------
    p%FRIMKP_PHI = FRIMKP_PHIb
    p%FRIMCP_PHI = FRIMCP_PHIb
    p%FRIMMP_PHI = FRIMMP_PHIb  
    ! AÑADO VERTICAL --------------------------------
    p%FRIMKP_V = FRIMKP_Vb
    p%FRIMCP_V = FRIMCP_Vb
    p%FRIMMP_V = FRIMMP_Vb     

     CALL AllocAry( p%PhiRbase, p%DOFL , p%DOFR - p%DOFI, 'p%PhiRbase', ErrStat2, ErrMsg2 ); if(Failed()) return
     p%PhiRbase = PhiRbase

    IF ( p%Nmodes > 0) THEN
       p%FMIMKP = FMIMKPb
       p%FMIMCP = FMIMCPb
       p%FMIMMP = FMIMMPb
    ! AÑADO GIRO --------------------------------
       p%FMIMKP_PHI = FMIMKP_PHIb
       p%FMIMCP_PHI = FMIMCP_PHIb
       p%FMIMMP_PHI = FMIMMP_PHIb   
     ! AÑADO VERTICAL --------------------------------
       p%FMIMKP_V = FMIMKP_Vb
       p%FMIMCP_V = FMIMCP_Vb
       p%FMIMMP_V = FMIMMP_Vb        
    END IF
   END IF

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SetParameters') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed
   
END SUBROUTINE SetParameters

!------------------------------------------------------------------------------------------------------

!> Allocate parameter arrays, based on the dimensions already set in the parameter data type.
SUBROUTINE AllocParameters(p, DOFM, ErrStat, ErrMsg)
   TYPE(SD_ParameterType), INTENT(INOUT)        :: p           ! Parameters
   INTEGER(IntKi), INTENT(  in)                 :: DOFM    
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2
   ! initialize error handling:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   ! for readability, we're going to keep track of the max ErrStat through SetErrStat() and not return until the end of this routine.
   
   CALL AllocAry( p%KBB,           TPdofL, TPdofL, 'p%KBB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%MBB,           TPdofL, TPdofL, 'p%MBB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%CBB,           TPdofL, TPdofL, 'p%CBB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')

   CALL AllocAry( p%TI,            p%DOFI,  6,     'p%TI',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%D1_14,         TPdofL, p%DOFL, 'p%D1_14',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%FY,            TPdofL,         'p%FY',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%PhiRb_TI,      p%DOFL, TPdofL, 'p%PhiRb_TI',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') 
   CALL AllocAry( p%FRIMKP,        TPdofL,         'p%FRIMKP',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')       
   CALL AllocAry( p%FRIMCP,        TPdofL,         'p%FRIMCP',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%FRIMMP,        TPdofL,         'p%FRIMMP',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') !AÑADO FM
   ! AÑADO GIRO
   CALL AllocAry( p%FRIMKP_PHI,        TPdofL,         'p%FRIMKP_PHI',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')       
   CALL AllocAry( p%FRIMCP_PHI,        TPdofL,         'p%FRIMCP_PHI',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%FRIMMP_PHI,        TPdofL,         'p%FRIMMP_PHI',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') 
   ! AÑADO VERTICAL
   CALL AllocAry( p%FRIMKP_V,        TPdofL,         'p%FRIMKP_V',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')       
   CALL AllocAry( p%FRIMCP_V,        TPdofL,         'p%FRIMCP_V',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%FRIMMP_V,        TPdofL,         'p%FRIMMP_V',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')    
 
if (p%Nmodes > 0 ) THEN  
   CALL AllocAry( p%MBM,           TPdofL, DOFM,   'p%MBM',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%MMB,           DOFM,   TPdofL, 'p%MMB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%CBM,           TPdofL, DOFM,   'p%CBM',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%CMB,           DOFM,   TPdofL, 'p%CMB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%CMM,           DOFM,     DOFM, 'p%CMM',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%NOmegaM2,      DOFM,           'p%NOmegaM2',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%N2OmegaMJDamp, DOFM,           'p%N2OmegaMJDamp', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%FX,            DOFM,           'p%FX',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%FMIMKP,        DOFM,           'p%FMIMKP',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')            
   CALL AllocAry( p%FMIMCP,        DOFM,           'p%FMIMCP',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') 
   CALL AllocAry( p%FMIMMP,        DOFM,           'p%FMIMMP',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')  !AÑADO FM   
   ! AÑADO GIRO
   CALL AllocAry( p%FMIMKP_PHI,        DOFM,           'p%FMIMKP_PHI',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')            
   CALL AllocAry( p%FMIMCP_PHI,        DOFM,           'p%FMIMCP_PHI',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') 
   CALL AllocAry( p%FMIMMP_PHI,        DOFM,           'p%FMIMMP_PHI',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')  
   ! AÑADO VERTICAL
   CALL AllocAry( p%FMIMKP_V,        DOFM,           'p%FMIMKP_V',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')            
   CALL AllocAry( p%FMIMCP_V,        DOFM,           'p%FMIMCP_V',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') 
   CALL AllocAry( p%FMIMMP_V,        DOFM,           'p%FMIMMP_V',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')     
   CALL AllocAry( p%C1_11,         TPdofL, DOFM,   'p%C1_11',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C1_12,         TPdofL, DOFM,   'p%C1_12',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%PhiM,          p%DOFL, DOFM,   'p%PhiM',          ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C2_61,         p%DOFL, DOFM,   'p%C2_61',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C2_62,         p%DOFL, DOFM,   'p%C2_62',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')    
   CALL AllocAry( p%D1_12,         TPdofL, TPdofL, 'p%D1_12',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')       
   CALL AllocAry( p%D1_13,         TPdofL, TPdofL, 'p%D1_13',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is p%MBB when p%NModes == 0 
   CALL AllocAry( p%D2_62,         p%DOFL, TPdofL, 'p%D2_62',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')       
   CALL AllocAry( p%D2_63,         p%DOFL, TPdofL, 'p%D2_63',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is p%PhiRb_TI when p%NModes == 0       
   CALL AllocAry( p%D2_64,         p%DOFL, p%DOFL, 'p%D2_64',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is zero when p%NModes == 0       
   CALL AllocAry( p%F2_61,         p%DOFL,         'p%F2_61',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is zero when p%NModes == 0
end if
                                   
   CALL AllocAry( p%IDI,           p%DOFI,               'p%IDI',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%IDR,           p%DOFR,               'p%IDR',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%IDL,           p%DOFL,               'p%IDL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%IDC,           p%DOFC,               'p%IDC',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%IDY,           p%DOFC+p%DOFI+p%DOFL, 'p%IDY',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
           
if ( p%SttcSolve ) THEN  
   CALL AllocAry( p%PhiL_T,        p%DOFL, p%DOFL, 'p%PhiL_T',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%PhiLInvOmgL2,  p%DOFL, p%DOFL, 'p%PhiLInvOmgL2',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%FGL,           p%DOFL,         'p%FGL',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')   
end if            
   
END SUBROUTINE AllocParameters

!------------------------------------------------------------------------------------------------------
!> Allocate parameter arrays, based on the dimensions already set in the parameter data type.
SUBROUTINE AllocMiscVars(p, Misc, ErrStat, ErrMsg)
   TYPE(SD_MiscVarType),    INTENT(INOUT)    :: Misc        ! Miscellaneous values, used to avoid local copies and/or multiple allocation/deallocation of same variables each call
   TYPE(SD_ParameterType),  INTENT(IN)       :: p           ! Parameters
   INTEGER(IntKi),          INTENT(  OUT)    :: ErrStat     ! Error status of the operation
   CHARACTER(*),            INTENT(  OUT)    :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   ! initialize error handling:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   ! for readability, we're going to keep track of the max ErrStat through SetErrStat() and not return until the end of this routine.
   CALL AllocAry( Misc%UFL,          p%DOFL,    'UFL',           ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UR_bar,       p%URbarL,  'UR_bar',        ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UR_bar_dot,   p%URbarL,  'UR_bar_dot',    ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UR_bar_dotdot,p%URbarL,  'UR_bar_dotdot', ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL,           p%DOFL,    'UL',            ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL_dot,       p%DOFL,    'UL_dot',        ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL_dotdot,    p%DOFL,    'UL_dotdot',     ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   
END SUBROUTINE AllocMiscVars

!------------------------------------------------------------------------------------------------------
!> Set the index arrays IDI, IDR, IDL, IDC, and IDY. 
SUBROUTINE SetIndexArrays(Init, p, ErrStat, ErrMsg)
   USE qsort_c_module, only: QsortC

   TYPE(SD_InitType),       INTENT(  IN)        :: Init        ! Input data for initialization routine
   TYPE(SD_ParameterType),  INTENT(INOUT)       :: p           ! Parameters   
   INTEGER(IntKi),          INTENT(  OUT)       :: ErrStat     ! Error status of the operation
   CHARACTER(*),            INTENT(  OUT)       :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER(IntKi)                               :: TempIDY(p%DOFC+p%DOFI+p%DOFL, 2)
   INTEGER(IntKi)                               :: IDT(Init%TDOF)
   INTEGER(IntKi)                               :: I, K  ! counters
   ErrStat = ErrID_None
   ErrMsg  = ""
         
   ! Index IDI for interface DOFs
   p%IDI = Init%IntFc(1:p%DOFI, 1)  !RRD interface DOFs
    
   ! Index IDC for constraint DOFs
   p%IDC = Init%BCs(1:p%DOFC, 1) !Constraint DOFs 
   
   ! Index IDR for IDR DOFs
   p%IDR(       1:p%DOFC ) = p%IDC  ! Constraint DOFs again
   p%IDR(p%DOFC+1:p%DOFR)  = p%IDI  ! IDR contains DOFs ofboundaries, constraints first then interface
   
   ! --- Index IDL for IDL DOFs
   ! first set the total DOFs:
   DO I = 1, Init%TDOF  !Total DOFs
      IDT(I) = I      
   ENDDO
   ! remove DOFs on the boundaries:
   DO I = 1, p%DOFR  !Boundary DOFs (Interface + Constraints)
      IDT(p%IDR(I)) = 0   !Set 0 wherever DOFs belong to boundaries
   ENDDO
   ! That leaves the internal DOFs:
   K = 0
   DO I = 1, Init%TDOF
      IF ( IDT(I) .NE. 0 ) THEN
         K = K+1
         p%IDL(K) = IDT(I)   !Internal DOFs
      ENDIF
   ENDDO   
   IF ( K /= p%DOFL ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = "SetIndexArrays: IDL or p%DOFL are the incorrect size."
      RETURN
   END IF
   
   ! --- Index IDY for all DOFs:
   ! set the second column of the temp array      
   DO I = 1, SIZE(TempIDY,1)
      TempIDY(I, 2) = I   ! this column will become the returned "key" (i.e., the original location in the array)
   ENDDO
   ! set the first column of the temp array      
   TempIDY(1:p%DOFI, 1) = p%IDI
   TempIDY(p%DOFI+1 : p%DOFI+p%DOFL, 1) = p%IDL
   TempIDY(p%DOFI+p%DOFL+1: p%DOFI+p%DOFL+p%DOFC, 1) = p%IDC
   ! sort based on the first column
   CALL QsortC( TempIDY )
   ! the second column is the key:
   p%IDY = TempIDY(:, 2)
   
END SUBROUTINE SetIndexArrays

!------------------------------------------------------------------------------------------------------
!>
SUBROUTINE Test_CB_Results(MBBt, MBMt, KBBt, OmegaM, DOFTP, DOFM, ErrStat, ErrMsg,Init,p)
   TYPE(SD_InitType),      INTENT(  in)                :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType), INTENT(inout)                :: p           ! Parameters
   INTEGER(IntKi)                                     :: DOFTP, DOFM
   REAL(ReKi)                                         :: MBBt(DOFTP, DOFTP)
   REAL(ReKi)                                         :: MBmt(DOFTP, DOFM)
   REAL(ReKi)                                         :: KBBt(DOFTP, DOFTP)
   REAL(ReKi)                                         :: OmegaM(DOFM)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER(IntKi) :: DOFT, NM, i
   REAL(ReKi), Allocatable     :: OmegaCB(:), PhiCB(:, :)
   REAL(ReKi), Allocatable     :: K(:, :)
   REAL(ReKi), Allocatable     :: M(:, :)
   Character(1024)             :: rootname
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   DOFT = DOFTP + DOFM
   NM = DOFT - 3
   Allocate( OmegaCB(NM), K(DOFT, DOFT), M(DOFT, DOFT), PhiCB(DOFT, NM) )
   K = 0.0
   M = 0.0
   OmegaCB = 0.0
   PhiCB = 0.0
   
   M(1:DOFTP, 1:DOFTP) = MBBt
   M(1:DOFTP, (DOFTP+1):DOFT ) = MBMt
   M((DOFTP+1):DOFT, 1:DOFTP ) = transpose(mbmt)

   DO i = 1, DOFM
      K(DOFTP+i, DOFTP+i) = OmegaM(i)*OmegaM(i)
      M(DOFTP+i, DOFTP+i) = 1.0
   ENDDO
      
   K(1:DOFTP, 1:DOFTP) = KBBt

   ! temporary rootname
   rootname = './test_assemble_C-B_out'
   
   CALL EigenSolve(K, M, DOFT, NM,.False.,Init,p, PhiCB, OmegaCB,  ErrStat, ErrMsg)
   IF ( ErrStat /= 0 ) RETURN  

END SUBROUTINE Test_CB_Results

!------------------------------------------------------------------------------------------------------
!> Take the input u LMesh and constructs the appropriate corresponding UFL vector
SUBROUTINE ConstructUFL( u, p, UFL )
   TYPE(SD_InputType),             INTENT(IN   )  :: u               ! Inputs
   TYPE(SD_ParameterType),         INTENT(IN   )  :: p               ! Parameters
   REAL(ReKi)                                     :: UFL(p%DOFL)
   INTEGER                                        :: I, J, StartDOF  ! integers for indexing into mesh and UFL

   ! note that p%DOFL = p%NNodes_L*6
   DO I = 1, p%NNodes_L   !Only interior nodes here     
      ! starting index in the master arrays for the current node    
      startDOF = (I-1)*6 + 1
      ! index into the Y2Mesh
      J  = p%NNodes_I + I
      ! Construct UFL array from the Force and Moment fields of the input mesh
      UFL ( startDOF   : startDOF + 2 ) = u%LMesh%Force (:,J)
      UFL ( startDOF+3 : startDOF + 5 ) = u%LMesh%Moment(:,J)
   END DO   

END SUBROUTINE

!------------------------------------------------------------------------------------------------------
!> Output the summary file    
SUBROUTINE OutSummary(Init, p, FEMparams,CBparams, ErrStat,ErrMsg)
   TYPE(SD_InitType),      INTENT(IN)     :: Init           ! Input data for initialization routine, this structure contains many variables needed for summary file
   TYPE(SD_ParameterType), INTENT(IN)     :: p              ! Parameters,this structure contains many variables needed for summary file
   TYPE(CB_MatArrays),     INTENT(IN)     :: CBparams       ! CB parameters that will be passed in for summary file use
   TYPE(FEM_MatArrays),    INTENT(IN)     :: FEMparams      ! FEM parameters that will be passed in for summary file use
   INTEGER(IntKi),         INTENT(OUT)    :: ErrStat        ! Error status of the operation
   CHARACTER(*),           INTENT(OUT)    :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   !LOCALS
   INTEGER(IntKi)         :: UnSum          ! unit number for this summary file
   INTEGER(IntKi)         :: ErrStat2       ! Temporary storage for local errors
   CHARACTER(ErrMsgLen)   :: ErrMsg2       ! Temporary storage for local errors
   CHARACTER(1024)        :: SummaryName    ! name of the SubDyn summary file
   INTEGER(IntKi)         :: i, j, k, propids(2)  !counter and temporary holders
   INTEGER(IntKi)         :: SDtoMeshIndx(Init%NNode)
   REAL(ReKi)             :: MRB(6,6)    !REDUCED SYSTEM Kmatrix, equivalent mass matrix
   REAL(ReKi)             :: XYZ1(3),XYZ2(3), DirCos(3,3), mlength !temporary arrays, member i-th direction cosine matrix (global to local) and member length
   CHARACTER(*),PARAMETER                 :: SectionDivide = '____________________________________________________________________________________________________'
   CHARACTER(*),PARAMETER                 :: SubSectionDivide = '__________'
   CHARACTER(2),  DIMENSION(6), PARAMETER :: MatHds= (/'X ', 'Y ', 'Z ', 'XX', 'YY', 'ZZ'/)  !Headers for the columns and rows of 6x6 matrices
   
   ErrStat = ErrID_None
   ErrMsg  = ""
    
   CALL SD_Y2Mesh_Mapping(p, SDtoMeshIndx )

   !-------------------------------------------------------------------------------------------------------------
   ! open txt file
   !-------------------------------------------------------------------------------------------------------------
   SummaryName = TRIM(Init%RootName)//'.sum'
   UnSum = -1            ! we haven't opened the summary file, yet.   

   CALL SDOut_OpenSum( UnSum, SummaryName, SD_ProgDesc, ErrStat2, ErrMsg2 )   
      CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Init' )
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnSum)
         RETURN
      END IF
      
   !-------------------------------------------------------------------------------------------------------------
   ! write discretized data to a txt file
   !-------------------------------------------------------------------------------------------------------------
!bjj: for debugging, i recommend using the p% versions of all these variables whenever possible in this summary file:
! (it helps in debugging)
   WRITE(UnSum, '(A)')  'Unless specified, units are consistent with Input units, [SI] system is advised.'
   WRITE(UnSum, '(A)') SectionDivide
      
   WRITE(UnSum, '()')    
   WRITE(UnSum, '(A,I6)')  'Number of nodes (NNodes):',Init%NNode
   WRITE(UnSum, '(A8,1x,A11,3(1x,A15))')  'Node No.', 'Y2Mesh Node',          'X (m)',           'Y (m)',           'Z (m)'         
   WRITE(UnSum, '(A8,1x,A11,3(1x,A15))')  '--------', '-----------', '---------------', '---------------', '---------------'
!   WRITE(UnSum, '(I8.0, E15.6,E15.6,E15.6)') (INT(Init%Nodes(i, 1)),(Init%Nodes(i, j), j = 2, JointsCol), i = 1, Init%NNode) !do not group the format or it won't work 3(E15.6) does not work !bjj???
   WRITE(UnSum, '('//Num2LStr(Init%NNode)//'(I8,3x,I9,'//Num2lstr(JointsCol-1)//'(1x,F15.4),:,/))') &
                          (NINT(Init%Nodes(i, 1)), SDtoMeshIndx(i), (Init%Nodes(i, j), j = 2, JointsCol), i = 1, Init%NNode)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'Number of elements (NElems):',Init%NElem
   WRITE(UnSum, '(A8,5(A10))')  'Elem No.',    'Node_I',     'Node_J',      'Prop_I',      'Prop_J',      'LPM_Flag'
   WRITE(UnSum, '(I8,I10,I10,I10,I10,I10)') ((p%Elems(i, j), j = 1, MembersCol), i = 1, Init%NElem)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'Number of Simplified LPM elements (NSLPMEl):',p%NSLPMEl
   
   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'Number of properties (NProps):',Init%NProp
   WRITE(UnSum, '(A8,5(A15))')  'Prop No.',     'YoungE',       'ShearG',       'MatDens',     'XsecD',      'XsecT'
   WRITE(UnSum, '(I8, E15.6,E15.6,E15.6,E15.6,E15.6 ) ') (NINT(Init%Props(i, 1)), (Init%Props(i, j), j = 2, 6), i = 1, Init%NProp)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'No. of Reaction DOFs:',p%NReact*6
   WRITE(UnSum, '(A, A6)')  'Reaction DOF_ID',      'LOCK'
   WRITE(UnSum, '(I10, I10)') ((Init%BCs(i, j), j = 1, 2), i = 1, p%NReact*6)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'No. of Interface DOFs:',p%DOFI
   WRITE(UnSum, '(A,A6)')  'Interface DOF ID',      'LOCK'
   WRITE(UnSum, '(I10, I10)') ((Init%IntFc(i, j), j = 1, 2), i = 1, p%DOFI)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'Number of concentrated masses (NCMass):',Init%NCMass
   WRITE(UnSum, '(A10,A15,A15,A15,A15)')  'JointCMass',     'Mass',         'JXX',             'JYY',             'JZZ'
   WRITE(UnSum, '(F10.0, E15.6,E15.6,E15.6,E15.6)') ((Init%Cmass(i, j), j = 1, 5), i = 1, Init%NCMass)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'Number of members',p%NMembers
   WRITE(UnSum, '(A,I6)')  'Number of nodes per beam member:', Init%Ndiv+1
   WRITE(UnSum, '(A9,A10,A10,A15,A16)')  'Member ID', 'Joint1_ID', 'Joint2_ID', 'Mass', 'Node IDs...'
   !WRITE(UnSum, '('//Num2LStr(Init%NDiv + 1 )//'(I6))') ((Init%MemberNodes(i, j), j = 1, Init%NDiv+1), i = 1, p%NMembers)
   DO i=1,p%NMembers
       !Calculate member mass here; this should really be done somewhere else, yet it is not used anywhere else
       !IT WILL HAVE TO BE MODIFIED FOR OTHER THAN CIRCULAR PIPE ELEMENTS
       propids=Init%Members(i,4:5)
       mlength=MemberLength(Init%Members(i,1),Init,ErrStat,ErrMsg)
       IF (ErrStat .EQ. ErrID_None) THEN
          IF (Init%Members(i,6) == 0) THEN ! beam
       WRITE(UnSum, '(I9,I10,I10, E15.6, A3,'//Num2LStr(Init%NDiv + 1 )//'(I6))')    Init%Members(i,1:3),                &
        MemberMass(Init%PropSets(propids(1),4),Init%PropSets(propids(1),5),Init%PropSets(propids(1),6),   &
                    Init%PropSets(propids(2),4),Init%PropSets(propids(2),5),Init%PropSets(propids(2),6), mlength, .TRUE.),  &
               ' ',(Init%MemberNodes(i, j), j = 1, Init%NDiv+1)
          ELSE ! SLPM
             WRITE(UnSum, '(I9,I10,I10,A3,'//Num2LStr(Init%NDiv + 1 )//'(I6))')    Init%Members(i,1:3),  &
                    ' ',(Init%MemberNodes(i, j), j = 1, Init%NDiv+1)
          ENDIF
       ELSE 
           RETURN
       ENDIF
   ENDDO   
   !-------------------------------------------------------------------------------------------------------------
   ! write Cosine matrix for all members to a txt file
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A, I6)') 'Direction Cosine Matrices for all Members: GLOBAL-2-LOCAL. No. of 3x3 matrices=', p%NMembers 
   WRITE(UnSum, '(A9,9(A15))')  'Member ID', 'DC(1,1)', 'DC(1,2)', 'DC(1,3)', 'DC(2,1)','DC(2,2)','DC(2,3)','DC(3,1)','DC(3,2)','DC(3,3)'
   DO i=1,p%NMembers
       !Find the right index in the Nodes array for the selected JointID. This is horrible, but I do not know how to implement this search in a more efficient way
       !The alternative would be to get an element that belongs to the member and use it with dircos
       
!BJJ:TODO:  DIDN'T we already calculate DirCos for each element? can't we use that here?       
       DO j=1,Init%NNode
           IF    ( NINT(Init%Nodes(j,1)) .EQ. Init%Members(i,2) )THEN 
                XYZ1=Init%Nodes(Init%Members(i,2),2:4)
           ELSEIF ( NINT(Init%Nodes(j,1)) .EQ. Init%Members(i,3) ) THEN 
                XYZ2=Init%Nodes(Init%Members(i,3),2:4)
           ENDIF
       ENDDO    
       CALL GetDirCos(XYZ1(1), XYZ1(2), XYZ1(3), XYZ2(1), XYZ2(2), XYZ2(3), DirCos, mlength, ErrStat, ErrMsg)
       DirCos=TRANSPOSE(DirCos) !This is now global to local
       WRITE(UnSum, '(I9,9(E15.6))') Init%Members(i,1), ((DirCos(k,j),j=1,3),k=1,3)
   ENDDO

   !-------------------------------------------------------------------------------------------------------------
   ! write Eigenvalues of full SYstem and CB reduced System
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'Eigenvalues'
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A, I6)') "FEM Eigenvalues [Hz]. Number of shown eigenvalues (total # of DOFs minus restrained nodes' DOFs):", FEMparams%NOmega 
   WRITE(UnSum, '(I6, e15.6)') ( i, FEMparams%Omega(i)/2.0/pi, i = 1, FEMparams%NOmega )

   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A, I6)') "CB Reduced Eigenvalues [Hz].  Number of retained modes' eigenvalues:", CBparams%DOFM 
   WRITE(UnSum, '(I6, e15.6)') ( i, CBparams%OmegaL(i)/2.0/pi, i = 1, CBparams%DOFM )  
    
   !-------------------------------------------------------------------------------------------------------------
   ! write Eigenvectors of full SYstem 
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A, I6)') ('FEM Eigenvectors ('//TRIM(Num2LStr(Init%TDOF))//' x '//TRIM(Num2LStr(FEMparams%NOmega))//&
                              ') [m or rad]. Number of shown eigenvectors (total # of DOFs minus restrained nodes'' DOFs):'), FEMparams%NOmega 
   WRITE(UnSum, '(6x,'//Num2LStr(FEMparams%NOmega)//'(I15))') (i, i = 1, FEMparams%NOmega  )!HEADERS
   WRITE(UnSum, '(I6,'//Num2LStr(FEMparams%NOmega)//'e15.6)') ( i, (FEMparams%Modes(i,j), j = 1, FEMparams%NOmega ),i = 1, Init%TDOF)

   IF (Init%CBMod) THEN  !with Craig-Bampton reduction
    
   !-------------------------------------------------------------------------------------------------------------
   ! write CB system matrices
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'CB Matrices (PhiM,PhiR) (no constraint applied)'
   
   WRITE(UnSum, '(A)') SubSectionDivide
   IF (CBparams%DOFM > 0) THEN
      CALL WrMatrix( CBparams%PhiL(:,1:CBparams%DOFM ), UnSum, 'e15.6', 'PhiM' ) 
   ELSE
      WRITE( UnSum, '(A,": ",A," x ",A)', IOSTAT=ErrStat ) "PhiM", TRIM(Num2LStr(p%DOFL)), '0' 
   END IF

   WRITE(UnSum, '(A)') SubSectionDivide
   CALL WrMatrix( CBparams%PhiR, UnSum, 'e15.6', 'PhiR' ) 
           
   !-------------------------------------------------------------------------------------------------------------
   ! write CB system KBBt and MBBt matrices, eq stiffness matrices of the entire substructure at the TP ref point
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') "SubDyn's Structure Equivalent Stiffness and Mass Matrices at the TP reference point (KBBt and MBBt)"
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A)') 'KBBt'  !Note p%KBB stores KBBt
   WRITE(UnSum, '(7(A15))') ' ', (MatHds(i), i = 1, 6   )
    !tried implicit loop unsuccessfully
    DO i=1,6
        WRITE(UnSum, '(A15, 6(e15.6))')   MatHds(i), (p%KBB(i,j), j = 1, 6)
    ENDDO    
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A)') ('MBBt')!Note p%MBB stores MBBt
   WRITE(UnSum, '(7(A15))') ' ', (MatHds(i), i = 1, 6   )
    DO i=1,6
        WRITE(UnSum, '(A15, 6(e15.6))')   MatHds(i), (p%MBB(i,j), j = 1, 6)
    ENDDO  
 
   MRB=matmul(TRANSPOSE(CBparams%TI2),matmul(CBparams%MBB,CBparams%TI2)) !Equivalent mass matrix of the rigid body
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'Rigid Body Equivalent Mass Matrix w.r.t. (0,0,0).'
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A)') 'MRB'
   WRITE(UnSum, '(7(A15))') ' ', (MatHds(i), i = 1, 6   )
   DO i=1,6
        WRITE(UnSum, '(A15, 6(e15.6))')   MatHds(i), (MRB(i,j), j = 1, 6)
   ENDDO 
   
   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,E15.6)')    "SubDyn's Total Mass (structural and non-structural)=", MRB(1,1) 
   WRITE(UnSum, '(A,3(E15.6))') "SubDyn's Total Mass CM coordinates (Xcm,Ycm,Zcm)   =", (/-MRB(3,5),-MRB(1,6), MRB(1,5)/) /MRB(1,1)        
   
!#ifdef SD_SUMMARY_DEBUG

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '**** Additional Debugging Information ****'

   ENDIF   !END, if Craig-Bampton reduction

   !-------------------------------------------------------------------------------------------------------------
   ! write assembed K C M to a txt file
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A, I6)') 'FULL FEM K and M matrices. TOTAL FEM TDOFs:', Init%TDOF 
   WRITE(UnSum, '(A)') ('Stiffness matrix K' )
   WRITE(UnSum, '(15x,'//TRIM(Num2LStr(Init%TDOF))//'(I15))')  (i, i = 1, Init%TDOF  )
   DO i=1,Init%TDOF
        WRITE(UnSum, '(I15, '//TRIM(Num2LStr(Init%TDOF))//'(e15.6))')   i, (Init%K(i, j), j = 1, Init%TDOF)
   ENDDO   

   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A)') ('Damping matrix C' )
   WRITE(UnSum, '(15x,'//TRIM(Num2LStr(Init%TDOF))//'(I15))')  (i, i = 1, Init%TDOF  )
   DO i=1,Init%TDOF
        WRITE(UnSum, '(I15, '//TRIM(Num2LStr(Init%TDOF))//'(e15.6))')   i, (Init%C(i, j), j = 1, Init%TDOF)
   ENDDO   
 
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A)') ('Mass matrix M' )
   WRITE(UnSum, '(15x,'//TRIM(Num2LStr(Init%TDOF))//'(I15))')  (i, i = 1, Init%TDOF  )
   DO i=1,Init%TDOF
        WRITE(UnSum, '(I15, '//TRIM(Num2LStr(Init%TDOF))//'(e15.6))')   i, (Init%M(i, j), j = 1, Init%TDOF)
   ENDDO  
   
   !-------------------------------------------------------------------------------------------------------------
   ! write assembed GRAVITY FORCE FG VECTOR.  gravity forces applied at each node of the full system
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'Gravity force vector FG applied at each node of the full system' 
   WRITE(UnSum, '(I6, e15.6)') (i, Init%FG(i), i = 1, Init%TDOF)

   IF (Init%CBMod) THEN  !with Craig-Bampton reduction
      
   !-------------------------------------------------------------------------------------------------------------
   ! write CB system matrices
   !-------------------------------------------------------------------------------------------------------------   
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'Additional CB Matrices (MBB,MBM,KBB) (no constraint applied)'
        
   WRITE(UnSum, '(A)') SubSectionDivide
   CALL WrMatrix( CBparams%MBB, UnSum, 'e15.6', 'MBB' ) 
    
   WRITE(UnSum, '(A)') SubSectionDivide
   IF ( CBparams%DOFM > 0 ) THEN
      CALL WrMatrix( CBparams%MBM, UnSum, 'e15.6', 'MBM' ) 
   ELSE
      WRITE( UnSum, '(A,": ",A," x ",A)', IOSTAT=ErrStat ) "MBM", '6', '0' 
   END IF
   
   WRITE(UnSum, '(A)') SubSectionDivide
   CALL WrMatrix( CBparams%KBB, UnSum, 'e15.6', 'KBB' ) 
    
   WRITE(UnSum, '(A)') SubSectionDivide
   CALL WrMatrix( CBparams%OmegaL**2, UnSum, 'e15.6','KMM (diagonal)' ) 

   ENDIF
   
   !-------------------------------------------------------------------------------------------------------------
   ! write TP TI matrix
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'TP refpoint Transformation Matrix TI '
   CALL WrMatrix( p%TI, UnSum, 'e15.6', 'TI' ) 
      
!#endif   
   
   CALL SDOut_CloseSum( UnSum, ErrStat, ErrMsg )  

END SUBROUTINE OutSummary

!------------------------------------------------------------------------------------------------------
!> This function calculates the length of a member 
FUNCTION MemberLength(MemberID,Init,ErrStat,ErrMsg)
    TYPE(SD_InitType), INTENT(IN)             :: Init         !< Input data for initialization routine, this structure contains many variables needed for summary file
    INTEGER(IntKi),    INTENT(IN)             :: MemberID     !< Member ID #
    REAL(ReKi)                                :: MemberLength !< Member Length
    INTEGER(IntKi),            INTENT(   OUT) :: ErrStat      !< Error status of the operation
    CHARACTER(*),              INTENT(   OUT) :: ErrMsg       !< Error message if ErrStat /= ErrID_None
    !Locals
    REAL(Reki)                    :: xyz1(3),xyz2(3)  ! Coordinates of joints in GLOBAL REF SYS
    INTEGER(IntKi)                :: i                ! Counter
    INTEGER(IntKi)                :: Joint1,Joint2    ! JointID
    CHARACTER(*), PARAMETER       :: RoutineName = 'MemberLength'
    ErrStat = ErrID_None
    ErrMsg  = ''
    MemberLength=0.0
    
    !Find the MemberID in the list
    DO i=1,SIZE(Init%Members, DIM=1)
        IF (Init%Members(i,1) .EQ. MemberID) THEN
           ! Find joints ID for this member
           Joint1 = FindNode(i,1); if (Joint1<0) return
           Joint2 = FindNode(i,2); if (Joint2<0) return
           xyz1= Init%Joints(Joint1,2:4)
           xyz2= Init%Joints(Joint2,2:4)
           MemberLength=SQRT( SUM((xyz2-xyz1)**2.) )
           if ( EqualRealNos(MemberLength, 0.0_ReKi) ) then 
               call SetErrStat(ErrID_Fatal,' Member with ID '//trim(Num2LStr(MemberID))//' has zero length!', ErrStat,ErrMsg,RoutineName);
               return
           endif
           return
       ENDIF
   ENDDO       
   call SetErrStat(ErrID_Fatal,' Member with ID '//trim(Num2LStr(MemberID))//' not found in member list!', ErrStat,ErrMsg,RoutineName);

contains
    !> Find JointID for node `iNode` (1 or 2) or member `iMember`
    integer(IntKi) function FindNode(iMember,iNode) result(j)
        integer(IntKi), intent(in) :: iMember !< Member index in Init%Members list
        integer(IntKi), intent(in) :: iNode   !< Node index, 1 or 2 for the member iMember
        logical  :: found
        found = .false.      
        j=1
        do while ( .not. found .and. j <= Init%NJoints )
            if (Init%Members(iMember, iNode+1) == nint(Init%Joints(j,1))) then ! Columns 2/3 for iNode 1/2
                found = .true.
                exit
            endif
            j = j + 1
        enddo 
        if (.not.found) then
            j=-1
            call SetErrStat(ErrID_Fatal,' Member '//trim(Num2LStr(iMember))//' has JointID'//trim(Num2LStr(iNode))//' = '//& 
                trim(Num2LStr(Init%Members(iMember,iNode+1)))//' which is not in the node list !', ErrStat,ErrMsg,RoutineName)
        endif
    end function

END FUNCTION MemberLength

!------------------------------------------------------------------------------------------------------
!> Calculate member mass, given properties at the ends, keep units consistent
!! For now it works only for circular pipes or for a linearly varying area
FUNCTION MemberMass(rho1,D1,t1,rho2,D2,t2,L,ctube)
    REAL(ReKi), INTENT(IN)                :: rho1,D1,t1,rho2,D2,t2 ,L       ! Density, OD and wall thickness for circular tube members at ends, Length of member
    !                                                     IF ctube=.FALSE. then D1/2=Area at end1/2, t1 and t2 are ignored
    REAL(ReKi)              :: MemberMass  !mass
    LOGICAL, INTENT(IN)                :: ctube          ! =TRUE for circular pipes, false elseshape
    !LOCALS
    REAL(ReKi)                ::a0,a1,a2,b0,b1,dd,dt  !temporary coefficients
    
    !Density allowed to vary linearly only
    b0=rho1
    b1=(rho2-rho1)/L
    !Here we will need to figure out what element it is for now circular pipes
        IF (ctube) THEN !circular tube
         a0=pi * (D1*t1-t1**2.)
         dt=t2-t1 !thickness variation
         dd=D2-D1 !OD variation
         a1=pi * ( dd*t1 + D1*dt -2.*t1*dt)/L 
         a2=pi * ( dd*dt-dt**2.)/L**2.
    
        ELSE  !linearly varying area
         a0=D1  !This is an area
         a1=(D2-D1)/L !Delta area
         a2=0.
    
        ENDIF
    MemberMass= b0*a0*L +(a0*b1+b0*a1)*L**2/2. + (b0*a2+b1*a1)*L**3/3 + a2*b1*L**4/4.!Integral of rho*A dz
      
END FUNCTION MemberMass

!------------------------------------------------------------------------------------------------------
!> Check whether MAT IS SYMMETRIC AND RETURNS THE MAXIMUM RELATIVE ERROR    
SUBROUTINE SymMatDebug(M,MAT)
    INTEGER(IntKi), INTENT(IN)                 :: M     ! Number of rows and columns
    REAL(ReKi),INTENT(IN)                      :: MAT(M ,M)    !matrix to be checked
    !LOCALS
    REAL(ReKi)                      :: Error,MaxErr    !element by element relative difference in (Transpose(MAT)-MAT)/MAT
    INTEGER(IntKi)                  ::  i, j, imax,jmax   !counter and temporary holders 

    MaxErr=0.
    imax=0
    jmax=0
    DO j=1,M
        DO i=1,M
            Error=MAT(i,j)-MAT(j,i)
            IF (MAT(i,j).NE.0) THEN
                Error=ABS(Error)/MAT(i,j)
            ENDIF    
            IF (Error > MaxErr) THEN
                imax=i
                jmax=j
                MaxErr=Error
            ENDIF    
        ENDDO
    ENDDO

   !--------------------------------------
   ! write discretized data to a txt file
   WRITE(*, '(A,e15.6)')  'Matrix Symmetry Check: Largest (abs) relative error is:', MaxErr
   WRITE(*, '(A,I4,I4)')  'Matrix Symmetry Check: (I,J)=', imax,jmax

END SUBROUTINE SymMatDebug

!-----------------------------------------------------------------------------------------------------
!> Saves in Ug and Ugdot the interpolated values of the input signal at current time t
!> Linear interporlation
SUBROUTINE InterpSeismicSignal(t, p, m)

      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables

      ! local variables

      REAL(ReKi)                                     :: frac
      INTEGER(IntKi)                                 :: I,Ifloor, Iceiling

      IFloor = FLOOR(t/p%SDDeltaTUg) + 1
      ICeiling = CEILING(t/p%SDDeltaTUg) + 1

      IF (IFloor == ICeiling) THEN

         m%Ug = p%UgData(Ifloor,2)
         m%Udotg = p%UgData(Ifloor,3)
         m%Uddotg = p%UgData(Ifloor,4)
         m%PHIg = p%UgData(Ifloor,5)
         m%PHIdotg = p%UgData(Ifloor,6)
         m%PHIddotg = p%UgData(Ifloor,7)
         m%Vg = p%UgData(Ifloor,8)
         m%Vdotg = p%UgData(Ifloor,9)
         m%Vddotg = p%UgData(Ifloor,10)

      ELSE

         I=Ifloor

         frac = ( t - p%UgData(I,1) ) / ( p%UgData(I+1,1) - p%UgData(I,1) )

         m%Ug     = p%UgData(I,2) + frac * ( p%UgData(I+1,2) - p%UgData(I,2) )
         m%Udotg  = p%UgData(I,3) + frac * ( p%UgData(I+1,3) - p%UgData(I,3) )
         m%Uddotg = p%UgData(I,4) + frac * ( p%UgData(I+1,4) - p%UgData(I,4) )
         m%PHIg     = p%UgData(I,5) + frac * ( p%UgData(I+1,5) - p%UgData(I,5) )
         m%PHIdotg  = p%UgData(I,6) + frac * ( p%UgData(I+1,6) - p%UgData(I,6) )
         m%PHIddotg = p%UgData(I,7) + frac * ( p%UgData(I+1,7) - p%UgData(I,7) )
         m%Vg     = p%UgData(I,8) + frac * ( p%UgData(I+1,8) - p%UgData(I,8) )
         m%Vdotg  = p%UgData(I,9) + frac * ( p%UgData(I+1,9) - p%UgData(I,9) )
         m%Vddotg = p%UgData(I,10) + frac * ( p%UgData(I+1,10) - p%UgData(I,10) )
         
         

      ENDIF

END SUBROUTINE InterpSeismicSignal

!! -------------------------------------------------------------!!
!! Rutines associated to FULL-FEM mode without modal reduction  !!
!! -------------------------------------------------------------!!

! Rutine needed to cast FULL FEM Matrices into state-space form later on
!------------------------------------------------------------------------------------------------------
SUBROUTINE FULLFEM(Init, p, FEMParams, ErrStat, ErrMsg)
   TYPE(SD_InitType),     INTENT(INOUT)      :: Init          ! Input data for initialization routine
   TYPE(SD_ParameterType),INTENT(INOUT)      :: p             ! Parameters
   TYPE(FEM_MatArrays)   ,INTENT(IN   )      :: FEMparams     ! FEM parameters
   INTEGER(IntKi),        INTENT(  OUT)      :: ErrStat       ! Error status of the operation
   CHARACTER(*),          INTENT(  OUT)      :: ErrMsg        ! Error message if ErrStat /= ErrID_None   

   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(ErrMsgLen)     :: ErrMsg2

   !Local variables
   REAL(ReKi), ALLOCATABLE  :: MII(:, :)
   REAL(ReKi), ALLOCATABLE  :: MLL(:, :)
   REAL(ReKi), ALLOCATABLE  :: MIL(:, :)
   REAL(ReKi), ALLOCATABLE  :: CII(:, :)
   REAL(ReKi), ALLOCATABLE  :: CLL(:, :)
   REAL(ReKi), ALLOCATABLE  :: CIL(:, :)
   REAL(ReKi), ALLOCATABLE  :: KII(:, :)
   REAL(ReKi), ALLOCATABLE  :: KLL(:, :)
   REAL(ReKi), ALLOCATABLE  :: KIL(:, :)
   REAL(ReKi), ALLOCATABLE  :: MIB(:, :)
   REAL(ReKi), ALLOCATABLE  :: MLB(:, :)

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! number of nodes:
   p%NNodes_I  = Init%NInterf                         ! Number of interface nodes
   p%NNodes_L  = Init%NNode - p%NReact - p%NNodes_I   ! Number of Interior nodes =(TDOF-DOFC-DOFI)/6 =  (6*Init%NNode - (p%NReact+p%NNodes_I)*6 ) / 6 = Init%NNode - p%NReact -p%NNodes_I

   !DOFS of interface
   !BJJ: TODO:  are these 6's actually TPdofL?   
   p%DOFI = p%NNodes_I*6
   p%DOFC = p%NReact*6
   p%DOFR = (p%NReact+p%NNodes_I)*6 ! = p%DOFC + p%DOFI
   p%DOFL = p%NNodes_L*6            ! = Init%TDOF - p%DOFR
         
   ! matrix dimension paramters
   p%URbarL = p%DOFI !=p%NNodes_I*6          ! Length of URbar array, subarray of Y2  : THIS MAY CHANGE IF SOME DOFS ARE NOT CONSTRAINED       
   
   ! Alloc matrices
   CALL AllocParametersFULLFEM(p, ErrStat2, ErrMsg2); 

   CALL AllocAry( MII,     p%DOFI, p%DOFI,        'matrix MII',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( MLL,     p%DOFL, p%DOFL,        'matrix MLL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( MIL,     p%DOFI, p%DOFL,        'matrix MIL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( CII,     p%DOFI, p%DOFI,        'matrix CII',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( CLL,     p%DOFL, p%DOFL,        'matrix CLL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( CIL,     p%DOFI, p%DOFL,        'matrix CIL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( KII,     p%DOFI, p%DOFI,        'matrix KII',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( KLL,     p%DOFL, p%DOFL,        'matrix KLL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( KIL,     p%DOFI, p%DOFL,        'matrix KIL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( MIB,     p%DOFI, p%DOFR-p%DOFI, 'matrix MIB',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( MLB,     p%DOFL, p%DOFR-p%DOFI, 'matrix MLB',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  

   ! After having performed Eigenvalue analysis, compute complete Damping Matrix with RAYLEIGH method. Caughey damping could be considered in the future
   ! If Lumped parameter elements are present, this damping matrix is added to the already existing one.  
   CALL RAYLEIGH(Init, p, FEMPARAMS)

   ! Set the index arrays p%IDI, p%IDR, p%IDL, p%IDC, and p%IDY. 
   CALL SetIndexArrays(Init, p, ErrStat2, ErrMsg2) ; if(Failed()) return

   ! Extract sub-matrices
   CALL BreakFEMSysMtrx(Init, p, MII, MIL, MLL, KII, KIL, KLL, CII, CIL, CLL, MIB, MLB)   
      
   ! Set p%TI
   CALL TrnsfTIFEM(Init, p%TI, p%DOFI, p%IDI, p%DOFR, p%IDR, ErrStat2, ErrMsg2); if(Failed()) return  

   !................................
   ! set values needed to calculate outputs and update states:
   !................................
   CALL SetParametersFULLFEM(Init, p, MII, MIL, MLL, KII, KIL, KLL, CII, CIL, CLL, MIB, MLB, ErrStat2, ErrMsg2)  
   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'Full FEM')
      
   CALL CleanUpFEM()

contains

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Full FEM') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUpFEM()
   end function Failed

   subroutine CleanUpFEM()
      IF(ALLOCATED(MII)  ) DEALLOCATE(MII) 
      IF(ALLOCATED(MLL)  ) DEALLOCATE(MLL) 
      IF(ALLOCATED(MIL)  ) DEALLOCATE(MIL)            
      IF(ALLOCATED(CII)  ) DEALLOCATE(CII) 
      IF(ALLOCATED(CLL)  ) DEALLOCATE(CLL) 
      IF(ALLOCATED(CIL)  ) DEALLOCATE(CIL)            
      IF(ALLOCATED(KII)  ) DEALLOCATE(KII) 
      IF(ALLOCATED(KLL)  ) DEALLOCATE(KLL) 
      IF(ALLOCATED(KIL)  ) DEALLOCATE(KIL)   
      IF(ALLOCATED(MIB)  ) DEALLOCATE(MIB) 
      IF(ALLOCATED(MLB)  ) DEALLOCATE(MLB)         
   end subroutine CleanUpFEM

END SUBROUTINE FULLFEM

!------------------------------------------------------------------------------------------------------
!>
SUBROUTINE BreakFEMSysMtrx(Init, p, MII, MIL, MLL, KII, KIL, KLL, CII, CIL, CLL, MIB, MLB)
   TYPE(SD_InitType),      INTENT(IN   )  :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType), INTENT(INOUT)  :: p  
   REAL(ReKi),             INTENT(  OUT)  :: MII(p%DOFI, p%DOFI)
   REAL(ReKi),             INTENT(  OUT)  :: MLL(p%DOFL, p%DOFL) 
   REAL(ReKi),             INTENT(  OUT)  :: MIL(p%DOFI, p%DOFL)
   REAL(ReKi),             INTENT(  OUT)  :: CII(p%DOFI, p%DOFI)
   REAL(ReKi),             INTENT(  OUT)  :: CLL(p%DOFL, p%DOFL) 
   REAL(ReKi),             INTENT(  OUT)  :: CIL(p%DOFI, p%DOFL)
   REAL(ReKi),             INTENT(  OUT)  :: KII(p%DOFI, p%DOFI)
   REAL(ReKi),             INTENT(  OUT)  :: KLL(p%DOFL, p%DOFL)
   REAL(ReKi),             INTENT(  OUT)  :: KIL(p%DOFI, p%DOFL)
   REAL(ReKi),             INTENT(  OUT)  :: MIB(p%DOFI, p%DOFR-p%DOFI)
   REAL(ReKi),             INTENT(  OUT)  :: MLB(p%DOFL, p%DOFR-p%DOFI) 

   ! local variables
   INTEGER(IntKi)          :: I, J, II, JJ
   
   DO I = 1, p%DOFI   
      II = p%IDI(I)
      p%FGI(I) = Init%FG(II)
      DO J = 1, p%DOFI
         JJ = p%IDI(J)
         MII(I, J) = Init%M(II, JJ)
         CII(I, J) = Init%C(II, JJ)
         KII(I, J) = Init%K(II, JJ)
      ENDDO
   ENDDO
   
   DO I = 1, p%DOFL
      II = p%IDL(I)
      p%FGL(I) = Init%FG(II)
      DO J = 1, p%DOFL
         JJ = p%IDL(J)
         MLL(I, J) = Init%M(II, JJ)
         CLL(I, J) = Init%C(II, JJ)
         KLL(I, J) = Init%K(II, JJ)
      ENDDO
   ENDDO
   
   DO I = 1, p%DOFI
      II = p%IDI(I)
      DO J = 1, p%DOFL
         JJ = p%IDL(J)
         MIL(I, J) = Init%M(II, JJ)
         CIL(I, J) = Init%C(II, JJ)
         KIL(I, J) = Init%K(II, JJ) 
      ENDDO                           
   ENDDO

   DO I = 1, p%DOFI
      II = p%IDI(I)
      DO J = 1, p%DOFR-p%DOFI
         JJ = p%IDC(J)
         MIB(I, J) = Init%M(II, JJ)
      ENDDO                           
   ENDDO

   DO I = 1, p%DOFL
      II = p%IDL(I)
      DO J = 1, p%DOFR-p%DOFI
         JJ = p%IDC(J)
         MLB(I, J) = Init%M(II, JJ)
      ENDDO                           
   ENDDO
      
END SUBROUTINE BreakFEMSysMtrx

!------------------------------------------------------------------------------------------------------
!>
SUBROUTINE TrnsfTIFEM(Init, TI, DOFI, IDI, DOFR, IDR, ErrStat, ErrMsg)
   TYPE(SD_InitType),      INTENT(IN   )  :: Init         ! Input data for initialization routine
   INTEGER(IntKi),         INTENT(IN   )  :: DOFI         ! # of DOFS of interface nodes
   INTEGER(IntKi),         INTENT(IN   )  :: DOFR         ! # of DOFS of restrained nodes (restraints and interface)
   INTEGER(IntKi),         INTENT(IN   )  :: IDI(DOFI)
   INTEGER(IntKi),         INTENT(IN   )  :: IDR(DOFR)
   REAL(ReKi),             INTENT(INOUT)  :: TI( DOFI,6)  ! matrix TI that relates the reduced matrix to the TP, 
   INTEGER(IntKi),         INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER                                :: I, di 
   INTEGER                                :: rmndr, n
   REAL(ReKi)                             :: dx, dy, dz
   
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   TI(:,:) = 0. !Initialize     
   DO I = 1, DOFI
      di = IDI(I)
      rmndr = MOD(di, 6)
      n = CEILING(di/6.0)
      
      dx = Init%Nodes(n, 2) - Init%TP_RefPoint(1)
      dy = Init%Nodes(n, 3) - Init%TP_RefPoint(2)
      dz = Init%Nodes(n, 4) - Init%TP_RefPoint(3)
      
      SELECT CASE (rmndr)
         CASE (1); TI(I, 1:6) = (/1.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi,       dz,      -dy/)
         CASE (2); TI(I, 1:6) = (/0.0_ReKi, 1.0_ReKi, 0.0_ReKi,      -dz, 0.0_ReKi,       dx/)
         CASE (3); TI(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi,       dy,      -dx, 0.0_ReKi/)
         CASE (4); TI(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi, 0.0_ReKi/)
         CASE (5); TI(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi, 0.0_ReKi/)
         CASE (0); TI(I, 1:6) = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, 1.0_ReKi/)
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMsg  = 'Error calculating transformation matrix TI '
            RETURN
         END SELECT
      
   ENDDO
   
END SUBROUTINE TrnsfTIFEM

!------------------------------------------------------------------------------------------------------

!> Allocate parameter arrays, based on the dimensions already set in the parameter data type.
SUBROUTINE AllocParametersFULLFEM(p, ErrStat, ErrMsg)
   TYPE(SD_ParameterType), INTENT(INOUT)        :: p           ! Parameters   
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2
   ! initialize error handling:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   CALL AllocAry( p%MII,          6,      6,        'matrix p%MII',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( p%MLL,     p%DOFL, p%DOFL,        'matrix p%MLL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( p%MIL,          6, p%DOFL,        'matrix p%MIL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%MLI,     p%DOFL,      6,        'matrix p%MLI',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%CII,          6,      6,        'matrix p%CII',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( p%CLL,     p%DOFL, p%DOFL,        'matrix p%CLL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( p%CIL,          6, p%DOFL,        'matrix p%CIL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%CLI,     p%DOFL,      6,        'matrix p%CLI',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%KII,          6,      6,        'matrix p%KII',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( p%KLL,     p%DOFL, p%DOFL,        'matrix p%KLL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( p%KIL,          6, p%DOFL,        'matrix p%KIL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%KLI,     p%DOFL,      6,        'matrix p%KIL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%MIB,     p%DOFI, p%DOFR-p%DOFI, 'matrix p%MIB',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')  
   CALL AllocAry( p%MLB,     p%DOFL, p%DOFR-p%DOFI, 'matrix p%MLB',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM') 

   CALL AllocAry( p%FSISKI,               6,        'p%FSISKI',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM') 
   CALL AllocAry( p%FSISCI,               6,        'p%FSISCI',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM') 
   CALL AllocAry( p%FSISMI,               6,        'p%FSISMI',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM') 
   CALL AllocAry( p%FSISKL,          p%DOFL,        'p%FSISKL',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM') 
   CALL AllocAry( p%FSISCL,          p%DOFL,        'p%FSISCL',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM') 
   CALL AllocAry( p%FSISML,          p%DOFL,        'p%FSISML',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM') 

   CALL AllocAry( p%FGL,             p%DOFL,        'p%FGL',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%FGI,             p%DOFI,        'p%FGI',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')

   CALL AllocAry( p%A_21,    p%DOFL, p%DOFL,        'p%A_21',          ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%A_22,    p%DOFL, p%DOFL,        'p%A_22',          ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%B_21,    p%DOFL,      6,        'p%B_21',          ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%B_22,    p%DOFL,      6,        'p%B_22',          ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%B_23,    p%DOFL,      6,        'p%B_23',          ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%B_24,    p%DOFL, p%DOFL,        'p%B_24',          ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%C1_11,        6, p%DOFL,        'p%C1_11',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%C1_12,        6, p%DOFL,        'p%C1_12',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%D1_11,        6,      6,        'p%D1_11',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%D1_12,        6,      6,        'p%D1_12',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%D1_13,        6,      6,        'p%D1_13',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%D1_14,        6, p%DOFL,        'p%D1_14',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%D1_15,        6, p%DOFI,        'p%D1_15',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL AllocAry( p%FX,              p%DOFL,        'p%FX',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM') 

   CALL AllocAry( p%TI,              p%DOFI,  6,     'p%TI',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')        
   CALL AllocAry( p%FY,              TPdofL,         'p%FY',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')    
                                   
   CALL AllocAry( p%IDI,           p%DOFI,               'p%IDI',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')        
   CALL AllocAry( p%IDR,           p%DOFR,               'p%IDR',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')        
   CALL AllocAry( p%IDL,           p%DOFL,               'p%IDL',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')        
   CALL AllocAry( p%IDC,           p%DOFC,               'p%IDC',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')        
   CALL AllocAry( p%IDY,           p%DOFC+p%DOFI+p%DOFL, 'p%IDY',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')        
           
END SUBROUTINE AllocParametersFULLFEM

!------------------------------------------------------------------------------------------------------
SUBROUTINE SetParametersFULLFEM(Init, p, MII, MIL, MLL, KII, KIL, KLL, CII, CIL, CLL, MIB, MLB,  ErrStat, ErrMsg)

   TYPE(SD_InitType),        INTENT(IN   )   :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType),   INTENT(INOUT)   :: p            ! Parameters
   REAL(ReKi),               INTENT(IN   )   :: MII(  p%DOFI, p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: MIL(  p%DOFI, p%DOFL)
   REAL(ReKi),               INTENT(IN   )   :: MLL(  p%DOFL, p%DOFL)
   REAL(ReKi),               INTENT(IN   )   :: KII(  p%DOFI, p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: KIL(  p%DOFI, p%DOFL)
   REAL(ReKi),               INTENT(IN   )   :: KLL(  p%DOFL, p%DOFL)
   REAL(ReKi),               INTENT(IN   )   :: CII(  p%DOFI, p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: CIL(  p%DOFI, p%DOFL)
   REAL(ReKi),               INTENT(IN   )   :: CLL(  p%DOFL, p%DOFL)
   REAL(ReKi),               INTENT(IN   )   :: MIB(  p%DOFI, p%DOFR - p%DOFI)
   REAL(ReKi),               INTENT(IN   )   :: MLB(  p%DOFL, p%DOFR - p%DOFI)


   INTEGER(IntKi),           INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   REAL(ReKi), ALLOCATABLE                   :: INVMLL(:, :)
   REAL(ReKi), ALLOCATABLE                   :: MIL_INVMLL(:, :)
   REAL(ReKi) , allocatable                  :: RI(:)             ! RI(p%DOFI) , influece vector on I dofs
   REAL(ReKi) , allocatable                  :: RL(:)             ! RL(p%DOFL) , influece vector on L dofs
   REAL(ReKi) , allocatable                  :: RB(:)             ! RL(p%DOFL) , influece vector on L dofs
   REAL(ReKi)                                :: TI_transpose(TPdofL,p%DOFI) !bjj: added this so we don't have to take the transpose 5+ times
   INTEGER(IntKi)                            :: I,J
   integer(IntKi)                            :: n                          ! size of jacobian in AM2 calculation
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   CHARACTER(*), PARAMETER                   :: RoutineName = 'SetParametersFULLFEM'
   
   ErrStat = ErrID_None 
   ErrMsg  = ''

   !Temporary computations that are used several times
      
   TI_transpose =  TRANSPOSE(p%TI) 

   CALL AllocAry( INVMLL, p%DOFL, p%DOFL, 'matrix INVMLL', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   INVMLL = INV(MLL,p%DOFL)
   
   !Matrices p%KLL, p%MII, etc., keep final tilde matrices in the FEM sistems of equations to be solved

   p%KLL = KLL   ! No transformation needed for II
   p%CLL = CLL
   p%MLL = MLL

   !Small matrices. Can be done efficiently with matmul
   p%MII = MATMUL( MATMUL( TI_transpose, MII ), p%TI) 
   p%KII = MATMUL( MATMUL( TI_transpose, KII ), p%TI) 
   p%CII = MATMUL( MATMUL( TI_transpose, CII ), p%TI) 

   CALL LAPACK_gemm( 'T', 'N', 1.0_ReKi, p%TI, MIL, 0.0_ReKi, p%MIL, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%MIL')
   CALL LAPACK_gemm( 'T', 'N', 1.0_ReKi, p%TI, KIL, 0.0_ReKi, p%KIL, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%KIL')
   CALL LAPACK_gemm( 'T', 'N', 1.0_ReKi, p%TI, CIL, 0.0_ReKi, p%CIL, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%CIL')

   p%MLI = TRANSPOSE( p%MIL )
   p%CLI = TRANSPOSE( p%CIL )
   p%KLI = TRANSPOSE( p%KIL )

   CALL AllocAry( MIL_INVMLL,  6, p%DOFL, 'matrix MIL_INVMLL', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FULLFEM')
   CALL LAPACK_gemm( 'N', 'N', 1.0_ReKi, p%MIL, INVMLL, 0.0_ReKi, MIL_INVMLL, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'MIL_INVMLL')

   ! --------------------------------------------------------------
   ! FORMULATION OF STATE-SPACE FORMULATION MATRICES AND VECTORS 
   ! -------------------------------------------------------------

   ! Matrix A

   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, INVMLL, p%KLL, 0.0_ReKi, p%A_21, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%A_21') 
   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, INVMLL, p%CLL, 0.0_ReKi, p%A_22, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%A_22') 

   ! Matrix B

   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, INVMLL, p%KLI, 0.0_ReKi, p%B_21, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%B_21') 
   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, INVMLL, p%CLI, 0.0_ReKi, p%B_22, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%B_22') 
   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, INVMLL, p%MLI, 0.0_ReKi, p%B_23, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%B_23')
   p%B_24 = INVMLL

   ! FX (except for the terms coming from the seismic input)

   p%FX = matmul(INVMLL , p%FGL)

   ! Matrix C1

   p%C1_11 = p%KIL
   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, MIL_INVMLL, p%KLL, 1.0_ReKi, p%C1_11, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%C1_11')  

   p%C1_12 = p%CIL
   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, MIL_INVMLL, p%CLL, 1.0_ReKi, p%C1_12, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%C1_12') 

   ! Matrix D1

   p%D1_11 = p%KII
   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, MIL_INVMLL, p%KLI, 1.0_ReKi, p%D1_11, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%D1_11') 

   p%D1_12 = p%CII
   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, MIL_INVMLL, p%CLI, 1.0_ReKi, p%D1_12, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%D1_12') 

   p%D1_13 = p%MII
   CALL LAPACK_gemm( 'N', 'N', -1.0_ReKi, MIL_INVMLL, p%MLI, 1.0_ReKi, p%D1_13, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%D1_13')

   p%D1_14 = MIL_INVMLL

   p%D1_15 = -TI_transpose

   ! FY (except for the terms coming from the seismic input)

   p%FY = - matmul(TI_transpose , p%FGI ) + matmul(MIL_INVMLL , p%FGL )

   ! Matrix C2 -> Neccesary items will be taken from A

   ! Matris D2 -> Neccesary items will be taken from B

   ! F2_61 -> Necessary item will be taken from FX      
                              
   !Now we should calculate a Jacobian used when AM2 is called and store in parameters    
    IF (p%IntMethod .EQ. 4) THEN       ! Allocate Jacobian if AM2 is requested & if there are states (p%qmL > 0)
       stop 'AM2 not implemented yet for FullFEM!' !The Jacobian should be reviewed before.
    END IF     

    IF (p%SeismicInp) THEN
      
    !................................
    ! set vectors related to seismic input (assumed horizontal only)
    ! The parts of the vectors that do not depend on t are set here
    !................................

    ! Set Influece vectors RR and RL
    
       CALL AllocAry( RL,  p%DOFL, 'Influence vector RL', ErrStat2, ErrMsg2); if(Failed()) return
       CALL AllocAry( RI,  p%DOFI, 'Influence vector RI', ErrStat2, ErrMsg2); if(Failed()) return
       CALL AllocAry( RB,  p%DOFR - p%DOFI, 'Influence vector RB', ErrStat2, ErrMsg2); if(Failed()) return
       CALL AllocAry( p%RRbase,  p%DOFR - p%DOFI, 'Influence vector RRbase', ErrStat2, ErrMsg2); if(Failed()) return

       RI = 0.0_ReKi
       RL = 0.0_ReKi
       RB = 0.0_ReKi
       J = 0

       DO I = 1,p%DOFI
          J = J + 1
          SELECT CASE (J)
             CASE (1); RI(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)  ! The argument of COS and SIN is radians. UgDir is in º
             CASE (2); RI(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
             CASE (6); J = 0
          END SELECT 
       END DO

       J = 0
       DO I = 1,p%DOFL
          J = J + 1
          SELECT CASE (J)
             CASE (1); RL(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)  ! The argument of COS and SIN is radians. UgDir is in º
             CASE (2); RL(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
             CASE (6); J = 0
          END SELECT 
       END DO

       J = 0
       DO I = 1,p%DOFR-P%DOFI
          J = J + 1
          SELECT CASE (J)
             CASE (1); RB(I) = COS((p%UgDir)*Pi_D/180.0_ReKi)  ! The argument of COS and SIN is radians. UgDir is in º
             CASE (2); RB(I) = SIN((p%UgDir)*Pi_D/180.0_ReKi)
             CASE (6); J = 0
          END SELECT 
       END DO

       p%RRbase = RB

       p%FSISKI = MATMUL(TI_transpose, ( MATMUL(KII , RI) + MATMUL(KIL, RL) ) )
       p%FSISCI = MATMUL(TI_transpose, ( MATMUL(CII , RI) + MATMUL(CIL, RL) ) )
       p%FSISMI = MATMUL(TI_transpose, ( MATMUL(MIB , RB) ) )

       p%FSISKL = MATMUL(transpose(KIL) , RI) + MATMUL(KLL, RL) 
       p%FSISCL = MATMUL(transpose(CIL) , RI) + MATMUL(CLL, RL)
       p%FSISML = MATMUL(MLB , RB) 

   
    DEALLOCATE(RI)
    DEALLOCATE(RL)
    DEALLOCATE(RB)
    END IF


    DEALLOCATE(INVMLL)
    DEALLOCATE(MIL_INVMLL)

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SetParametersFULLFEM') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed
   
END SUBROUTINE SetParametersFULLFEM


! -------------------------------------------------------------------------------
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A,n) result(Ainv)
  integer(IntKi)            , intent(in) :: n
  real(ReKi), dimension(:,:), intent(in) :: A
  real(ReKi), dimension(n,n) :: Ainv

  real(ReKi), dimension(n) :: work  ! work array for LAPACK
  integer, dimension(n) :: ipiv   ! pivot indices
  integer :: info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv

End Module SubDyn
