!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qes_write_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes-1.0
  !
  USE qes_types_module
  USE FoX_wxml
  !
  IMPLICIT NONE
  !
  PUBLIC
  !
  INTERFACE qes_write
    MODULE PROCEDURE qes_write_general_info
    MODULE PROCEDURE qes_write_parallel_info
    MODULE PROCEDURE qes_write_input
    MODULE PROCEDURE qes_write_step
    MODULE PROCEDURE qes_write_output
    MODULE PROCEDURE qes_write_control_variables
    MODULE PROCEDURE qes_write_xml_format
    MODULE PROCEDURE qes_write_creator
    MODULE PROCEDURE qes_write_created
    MODULE PROCEDURE qes_write_atomic_species
    MODULE PROCEDURE qes_write_species
    MODULE PROCEDURE qes_write_atomic_structure
    MODULE PROCEDURE qes_write_atomic_positions
    MODULE PROCEDURE qes_write_atom
    MODULE PROCEDURE qes_write_wyckoff_positions
    MODULE PROCEDURE qes_write_cell
    MODULE PROCEDURE qes_write_dft
    MODULE PROCEDURE qes_write_hybrid
    MODULE PROCEDURE qes_write_qpoint_grid
    MODULE PROCEDURE qes_write_dftU
    MODULE PROCEDURE qes_write_HubbardCommon
    MODULE PROCEDURE qes_write_HubbardJ
    MODULE PROCEDURE qes_write_starting_ns
    MODULE PROCEDURE qes_write_Hubbard_ns
    MODULE PROCEDURE qes_write_vdW
    MODULE PROCEDURE qes_write_spin
    MODULE PROCEDURE qes_write_bands
    MODULE PROCEDURE qes_write_smearing
    MODULE PROCEDURE qes_write_occupations
    MODULE PROCEDURE qes_write_basis
    MODULE PROCEDURE qes_write_basis_set
    MODULE PROCEDURE qes_write_basisSetItem
    MODULE PROCEDURE qes_write_reciprocal_lattice
    MODULE PROCEDURE qes_write_electron_control
    MODULE PROCEDURE qes_write_k_points_IBZ
    MODULE PROCEDURE qes_write_monkhorst_pack
    MODULE PROCEDURE qes_write_k_point
    MODULE PROCEDURE qes_write_ion_control
    MODULE PROCEDURE qes_write_bfgs
    MODULE PROCEDURE qes_write_md
    MODULE PROCEDURE qes_write_cell_control
    MODULE PROCEDURE qes_write_symmetry_flags
    MODULE PROCEDURE qes_write_boundary_conditions
    MODULE PROCEDURE qes_write_esm
    MODULE PROCEDURE qes_write_ekin_functional
    MODULE PROCEDURE qes_write_spin_constraints
    MODULE PROCEDURE qes_write_electric_field
    MODULE PROCEDURE qes_write_gate_settings
    MODULE PROCEDURE qes_write_atomic_constraints
    MODULE PROCEDURE qes_write_atomic_constraint
    MODULE PROCEDURE qes_write_inputOccupations
    MODULE PROCEDURE qes_write_outputElectricField
    MODULE PROCEDURE qes_write_BerryPhaseOutput
    MODULE PROCEDURE qes_write_dipoleOutput
    MODULE PROCEDURE qes_write_finiteFieldOut
    MODULE PROCEDURE qes_write_polarization
    MODULE PROCEDURE qes_write_ionicPolarization
    MODULE PROCEDURE qes_write_electronicPolarization
    MODULE PROCEDURE qes_write_phase
    MODULE PROCEDURE qes_write_gateInfo
    MODULE PROCEDURE qes_write_convergence_info
    MODULE PROCEDURE qes_write_scf_conv
    MODULE PROCEDURE qes_write_opt_conv
    MODULE PROCEDURE qes_write_algorithmic_info
    MODULE PROCEDURE qes_write_symmetries
    MODULE PROCEDURE qes_write_symmetry
    MODULE PROCEDURE qes_write_equivalent_atoms
    MODULE PROCEDURE qes_write_info
    MODULE PROCEDURE qes_write_outputPBC
    MODULE PROCEDURE qes_write_magnetization
    MODULE PROCEDURE qes_write_total_energy
    MODULE PROCEDURE qes_write_band_structure
    MODULE PROCEDURE qes_write_ks_energies
    MODULE PROCEDURE qes_write_closed
    MODULE PROCEDURE qes_write_vector
    MODULE PROCEDURE qes_write_integerVector
    MODULE PROCEDURE qes_write_matrix
    MODULE PROCEDURE qes_write_integerMatrix
    MODULE PROCEDURE qes_write_scalarQuantity
  END INTERFACE qes_write
  !
  CONTAINS
  !
  
   SUBROUTINE qes_write_general_info(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(general_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_xml_format (xp, obj%xml_format)
     CALL qes_write_creator (xp, obj%creator)
     CALL qes_write_created (xp, obj%created)
     CALL xml_addCharacters(xp, obj%job)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_general_info

   SUBROUTINE qes_write_parallel_info(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(parallel_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%nprocs)
     CALL xml_addCharacters(xp, obj%nthreads)
     CALL xml_addCharacters(xp, obj%ntasks)
     CALL xml_addCharacters(xp, obj%nbgrp)
     CALL xml_addCharacters(xp, obj%npool)
     CALL xml_addCharacters(xp, obj%ndiag)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_parallel_info

   SUBROUTINE qes_write_input(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(input_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_control_variables (xp, obj%control_variables)
     CALL qes_write_atomic_species (xp, obj%atomic_species)
     CALL qes_write_atomic_structure (xp, obj%atomic_structure)
     CALL qes_write_dft (xp, obj%dft)
     CALL qes_write_spin (xp, obj%spin)
     CALL qes_write_bands (xp, obj%bands)
     CALL qes_write_basis (xp, obj%basis)
     CALL qes_write_electron_control (xp, obj%electron_control)
     CALL qes_write_k_points_IBZ (xp, obj%k_points_IBZ)
     CALL qes_write_ion_control (xp, obj%ion_control)
     CALL qes_write_cell_control (xp, obj%cell_control)
     IF (obj%symmetry_flags_ispresent) THEN
        CALL qes_write_symmetry_flags (xp, obj%symmetry_flags)
     END IF
     IF (obj%boundary_conditions_ispresent) THEN
        CALL qes_write_boundary_conditions (xp, obj%boundary_conditions)
     END IF
     IF (obj%ekin_functional_ispresent) THEN
        CALL qes_write_ekin_functional (xp, obj%ekin_functional)
     END IF
     IF (obj%external_atomic_forces_ispresent) THEN
        CALL qes_write_matrix (xp, obj%external_atomic_forces)
     END IF
     IF (obj%free_positions_ispresent) THEN
        CALL qes_write_integerMatrix (xp, obj%free_positions)
     END IF
     IF (obj%starting_atomic_velocities_ispresent) THEN
        CALL qes_write_matrix (xp, obj%starting_atomic_velocities)
     END IF
     IF (obj%electric_field_ispresent) THEN
        CALL qes_write_electric_field (xp, obj%electric_field)
     END IF
     IF (obj%atomic_constraints_ispresent) THEN
        CALL qes_write_atomic_constraints (xp, obj%atomic_constraints)
     END IF
     IF (obj%spin_constraints_ispresent) THEN
        CALL qes_write_spin_constraints (xp, obj%spin_constraints)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_input

   SUBROUTINE qes_write_step(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(step_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'n_step', obj%n_step )
     CALL qes_write_scf_conv (xp, obj%scf_conv)
     CALL qes_write_atomic_structure (xp, obj%atomic_structure)
     CALL qes_write_total_energy (xp, obj%total_energy)
     CALL qes_write_matrix (xp, obj%forces)
     IF (obj%stress_ispresent) THEN
        CALL qes_write_matrix (xp, obj%stress)
     END IF
     IF (obj%FCP_force_ispresent) THEN
        CALL xml_addCharacters(xp, obj%FCP_force)
     END IF
     IF (obj%FCP_tot_charge_ispresent) THEN
        CALL xml_addCharacters(xp, obj%FCP_tot_charge)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_step

   SUBROUTINE qes_write_output(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(output_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_convergence_info (xp, obj%convergence_info)
     CALL qes_write_algorithmic_info (xp, obj%algorithmic_info)
     CALL qes_write_atomic_species (xp, obj%atomic_species)
     CALL qes_write_atomic_structure (xp, obj%atomic_structure)
     IF (obj%symmetries_ispresent) THEN
        CALL qes_write_symmetries (xp, obj%symmetries)
     END IF
     CALL qes_write_basis_set (xp, obj%basis_set)
     CALL qes_write_dft (xp, obj%dft)
     IF (obj%boundary_conditions_ispresent) THEN
        CALL qes_write_outputPBC (xp, obj%boundary_conditions)
     END IF
     CALL qes_write_magnetization (xp, obj%magnetization)
     CALL qes_write_total_energy (xp, obj%total_energy)
     CALL qes_write_band_structure (xp, obj%band_structure)
     IF (obj%forces_ispresent) THEN
        CALL qes_write_matrix (xp, obj%forces)
     END IF
     IF (obj%stress_ispresent) THEN
        CALL qes_write_matrix (xp, obj%stress)
     END IF
     IF (obj%electric_field_ispresent) THEN
        CALL qes_write_outputElectricField (xp, obj%electric_field)
     END IF
     IF (obj%FCP_force_ispresent) THEN
        CALL xml_addCharacters(xp, obj%FCP_force)
     END IF
     IF (obj%FCP_tot_charge_ispresent) THEN
        CALL xml_addCharacters(xp, obj%FCP_tot_charge)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_output

   SUBROUTINE qes_write_control_variables(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(control_variables_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%title)
     CALL xml_addCharacters(xp, obj%calculation)
     CALL xml_addCharacters(xp, obj%restart_mode)
     CALL xml_addCharacters(xp, obj%prefix)
     CALL xml_addCharacters(xp, obj%pseudo_dir)
     CALL xml_addCharacters(xp, obj%outdir)
     CALL xml_addCharacters(xp, obj%stress)
     CALL xml_addCharacters(xp, obj%forces)
     CALL xml_addCharacters(xp, obj%wf_collect)
     CALL xml_addCharacters(xp, obj%disk_io)
     CALL xml_addCharacters(xp, obj%max_seconds)
     IF (obj%nstep_ispresent) THEN
        CALL xml_addCharacters(xp, obj%nstep)
     END IF
     CALL xml_addCharacters(xp, obj%etot_conv_thr)
     CALL xml_addCharacters(xp, obj%forc_conv_thr)
     CALL xml_addCharacters(xp, obj%press_conv_thr)
     CALL xml_addCharacters(xp, obj%verbosity)
     CALL xml_addCharacters(xp, obj%print_every)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_control_variables

   SUBROUTINE qes_write_xml_format(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(xml_format_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'NAME', TRIM(obj%NAME) )
     CALL xml_addAttribute(xp, 'VERSION', TRIM(obj%VERSION) )
     CALL xml_AddCharacters(xp, obj%xml_format)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_xml_format

   SUBROUTINE qes_write_creator(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(creator_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'NAME', TRIM(obj%NAME) )
     CALL xml_addAttribute(xp, 'VERSION', TRIM(obj%VERSION) )
     CALL xml_AddCharacters(xp, obj%creator)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_creator

   SUBROUTINE qes_write_created(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(created_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'DATE', TRIM(obj%DATE) )
     CALL xml_addAttribute(xp, 'TIME', TRIM(obj%TIME) )
     CALL xml_AddCharacters(xp, obj%created)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_created

   SUBROUTINE qes_write_atomic_species(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_species_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'ntyp', obj%ntyp )
     IF (obj%pseudo_dir_ispresent) CALL xml_addAttribute(xp, 'pseudo_dir', TRIM(obj%pseudo_dir) )
     DO i = 1, obj%ndim_species
        CALL qes_write_species(xp, obj%species(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_species

   SUBROUTINE qes_write_species(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(species_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'name', TRIM(obj%name) )
     IF (obj%mass_ispresent) THEN
        CALL xml_addCharacters(xp, obj%mass)
     END IF
     CALL xml_addCharacters(xp, obj%pseudo_file)
     IF (obj%starting_magnetization_ispresent) THEN
        CALL xml_addCharacters(xp, obj%starting_magnetization)
     END IF
     IF (obj%spin_teta_ispresent) THEN
        CALL xml_addCharacters(xp, obj%spin_teta)
     END IF
     IF (obj%spin_phi_ispresent) THEN
        CALL xml_addCharacters(xp, obj%spin_phi)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_species

   SUBROUTINE qes_write_atomic_structure(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_structure_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'nat', obj%nat )
     IF (obj%alat_ispresent) CALL xml_addAttribute(xp, 'alat', obj%alat )
     IF (obj%bravais_index_ispresent) CALL xml_addAttribute(xp, 'bravais_index', obj%bravais_index )
     IF (obj%atomic_positions_ispresent) THEN
        CALL qes_write_atomic_positions (xp, obj%atomic_positions)
     END IF
     IF (obj%wyckoff_positions_ispresent) THEN
        CALL qes_write_wyckoff_positions (xp, obj%wyckoff_positions)
     END IF
     IF (obj%crystal_positions_ispresent) THEN
        CALL qes_write_atomic_positions (xp, obj%crystal_positions)
     END IF
     CALL qes_write_cell (xp, obj%cell)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_structure

   SUBROUTINE qes_write_atomic_positions(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_positions_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     DO i = 1, obj%ndim_atom
        CALL qes_write_atom(xp, obj%atom(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_positions

   SUBROUTINE qes_write_atom(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atom_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'name', TRIM(obj%name) )
     IF (obj%position_ispresent) CALL xml_addAttribute(xp, 'position', TRIM(obj%position) )
     IF (obj%index_ispresent) CALL xml_addAttribute(xp, 'index', obj%index )
     CALL xml_AddCharacters(xp, obj%atom)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atom

   SUBROUTINE qes_write_wyckoff_positions(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(wyckoff_positions_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'space_group', obj%space_group )
     IF (obj%more_options_ispresent) CALL xml_addAttribute(xp, 'more_options', TRIM(obj%more_options) )
     DO i = 1, obj%ndim_atom
        CALL qes_write_atom(xp, obj%atom(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_wyckoff_positions

   SUBROUTINE qes_write_cell(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cell_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%a1)
     CALL xml_addCharacters(xp, obj%a2)
     CALL xml_addCharacters(xp, obj%a3)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cell

   SUBROUTINE qes_write_dft(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(dft_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%functional)
     IF (obj%hybrid_ispresent) THEN
        CALL qes_write_hybrid (xp, obj%hybrid)
     END IF
     IF (obj%dftU_ispresent) THEN
        CALL qes_write_dftU (xp, obj%dftU)
     END IF
     IF (obj%vdW_ispresent) THEN
        CALL qes_write_vdW (xp, obj%vdW)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_dft

   SUBROUTINE qes_write_hybrid(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(hybrid_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_qpoint_grid (xp, obj%qpoint_grid)
     CALL xml_addCharacters(xp, obj%ecutfock)
     CALL xml_addCharacters(xp, obj%exx_fraction)
     CALL xml_addCharacters(xp, obj%screening_parameter)
     CALL xml_addCharacters(xp, obj%exxdiv_treatment)
     CALL xml_addCharacters(xp, obj%x_gamma_extrapolation)
     CALL xml_addCharacters(xp, obj%ecutvcut)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_hybrid

   SUBROUTINE qes_write_qpoint_grid(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(qpoint_grid_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'nqx1', obj%nqx1 )
     CALL xml_addAttribute(xp, 'nqx2', obj%nqx2 )
     CALL xml_addAttribute(xp, 'nqx3', obj%nqx3 )
     CALL xml_AddCharacters(xp, obj%qpoint_grid)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_qpoint_grid

   SUBROUTINE qes_write_dftU(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(dftU_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%lda_plus_u_kind_ispresent) THEN
        CALL xml_addCharacters(xp, obj%lda_plus_u_kind)
     END IF
     IF (obj%Hubbard_U_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_U
           CALL qes_write_HubbardCommon(xp, obj%Hubbard_U(i) )
        END DO
     END IF
     IF (obj%Hubbard_J0_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_J0
           CALL qes_write_HubbardCommon(xp, obj%Hubbard_J0(i) )
        END DO
     END IF
     IF (obj%Hubbard_alpha_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_alpha
           CALL qes_write_HubbardCommon(xp, obj%Hubbard_alpha(i) )
        END DO
     END IF
     IF (obj%Hubbard_beta_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_beta
           CALL qes_write_HubbardCommon(xp, obj%Hubbard_beta(i) )
        END DO
     END IF
     IF (obj%Hubbard_J_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_J
           CALL qes_write_HubbardJ(xp, obj%Hubbard_J(i) )
        END DO
     END IF
     IF (obj%starting_ns_ispresent) THEN
        DO i = 1, obj%ndim_starting_ns
           CALL qes_write_starting_ns(xp, obj%starting_ns(i) )
        END DO
     END IF
     IF (obj%Hubbard_ns_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_ns
           CALL qes_write_Hubbard_ns(xp, obj%Hubbard_ns(i) )
        END DO
     END IF
     IF (obj%U_projection_type_ispresent) THEN
        CALL xml_addCharacters(xp, obj%U_projection_type)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_dftU

   SUBROUTINE qes_write_HubbardCommon(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(HubbardCommon_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     CALL xml_AddCharacters(xp, obj%HubbardCommon)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_HubbardCommon

   SUBROUTINE qes_write_HubbardJ(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(HubbardJ_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     CALL xml_AddCharacters(xp, obj%HubbardJ)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_HubbardJ

   SUBROUTINE qes_write_starting_ns(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(starting_ns_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     CALL xml_addAttribute(xp, 'spin', obj%spin )
     CALL xml_addAttribute(xp, 'size', obj%size )
     CALL xml_AddCharacters(xp, obj%starting_ns, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_starting_ns

   SUBROUTINE qes_write_Hubbard_ns(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(Hubbard_ns_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     CALL xml_addAttribute(xp, 'spin', obj%spin )
     CALL xml_addAttribute(xp, 'index', obj%index )
     CALL xml_addAttribute(xp, 'rank', obj%rank )
     CALL xml_addAttribute(xp, 'dims', obj%dims )
     CALL xml_addAttribute(xp, 'order', TRIM(obj%order) )
       CALL xml_addNewLine(xp)
        DO i = 1, obj%dims(2)
           CALL xml_AddCharacters(xp, obj%Hubbard_ns((i-1)*obj%dims(1)+1: i*obj%dims(1)), fmt ='s16')
           CALL xml_addNewLine(xp)
        END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_Hubbard_ns

   SUBROUTINE qes_write_vdW(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(vdW_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%vdw_corr)
     IF (obj%non_local_term_ispresent) THEN
        CALL xml_addCharacters(xp, obj%non_local_term)
     END IF
     IF (obj%london_s6_ispresent) THEN
        CALL xml_addCharacters(xp, obj%london_s6)
     END IF
     IF (obj%ts_vdw_econv_thr_ispresent) THEN
        CALL xml_addCharacters(xp, obj%ts_vdw_econv_thr)
     END IF
     IF (obj%ts_vdw_isolated_ispresent) THEN
        CALL xml_addCharacters(xp, obj%ts_vdw_isolated)
     END IF
     IF (obj%london_rcut_ispresent) THEN
        CALL xml_addCharacters(xp, obj%london_rcut)
     END IF
     IF (obj%xdm_a1_ispresent) THEN
        CALL xml_addCharacters(xp, obj%xdm_a1)
     END IF
     IF (obj%xdm_a2_ispresent) THEN
        CALL xml_addCharacters(xp, obj%xdm_a2)
     END IF
     IF (obj%london_c6_ispresent) THEN
        DO i = 1, obj%ndim_london_c6
           CALL qes_write_HubbardCommon(xp, obj%london_c6(i) )
        END DO
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_vdW

   SUBROUTINE qes_write_spin(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(spin_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%lsda)
     CALL xml_addCharacters(xp, obj%noncolin)
     CALL xml_addCharacters(xp, obj%spinorbit)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_spin

   SUBROUTINE qes_write_bands(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(bands_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%nbnd_ispresent) THEN
        CALL xml_addCharacters(xp, obj%nbnd)
     END IF
     IF (obj%smearing_ispresent) THEN
        CALL qes_write_smearing (xp, obj%smearing)
     END IF
     IF (obj%tot_charge_ispresent) THEN
        CALL xml_addCharacters(xp, obj%tot_charge)
     END IF
     IF (obj%tot_magnetization_ispresent) THEN
        CALL xml_addCharacters(xp, obj%tot_magnetization)
     END IF
     CALL qes_write_occupations (xp, obj%occupations)
     IF (obj%inputOccupations_ispresent) THEN
        DO i = 1, obj%ndim_inputOccupations
           CALL qes_write_inputOccupations(xp, obj%inputOccupations(i) )
        END DO
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_bands

   SUBROUTINE qes_write_smearing(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(smearing_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'degauss', obj%degauss )
     CALL xml_AddCharacters(xp, obj%smearing)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_smearing

   SUBROUTINE qes_write_occupations(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(occupations_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%spin_ispresent) CALL xml_addAttribute(xp, 'spin', obj%spin )
     CALL xml_AddCharacters(xp, obj%occupations)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_occupations

   SUBROUTINE qes_write_basis(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(basis_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%gamma_only_ispresent) THEN
        CALL xml_addCharacters(xp, obj%gamma_only)
     END IF
     CALL xml_addCharacters(xp, obj%ecutwfc)
     IF (obj%ecutrho_ispresent) THEN
        CALL xml_addCharacters(xp, obj%ecutrho)
     END IF
     IF (obj%fft_grid_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_grid)
     END IF
     IF (obj%fft_smooth_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_smooth)
     END IF
     IF (obj%fft_box_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_box)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_basis

   SUBROUTINE qes_write_basis_set(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(basis_set_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%gamma_only_ispresent) THEN
        CALL xml_addCharacters(xp, obj%gamma_only)
     END IF
     CALL xml_addCharacters(xp, obj%ecutwfc)
     IF (obj%ecutrho_ispresent) THEN
        CALL xml_addCharacters(xp, obj%ecutrho)
     END IF
     CALL qes_write_basisSetItem (xp, obj%fft_grid)
     IF (obj%fft_smooth_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_smooth)
     END IF
     IF (obj%fft_box_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_box)
     END IF
     CALL xml_addCharacters(xp, obj%ngm)
     IF (obj%ngms_ispresent) THEN
        CALL xml_addCharacters(xp, obj%ngms)
     END IF
     CALL xml_addCharacters(xp, obj%npwx)
     CALL qes_write_reciprocal_lattice (xp, obj%reciprocal_lattice)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_basis_set

   SUBROUTINE qes_write_basisSetItem(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(basisSetItem_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'nr1', obj%nr1 )
     CALL xml_addAttribute(xp, 'nr2', obj%nr2 )
     CALL xml_addAttribute(xp, 'nr3', obj%nr3 )
     CALL xml_AddCharacters(xp, obj%basisSetItem)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_basisSetItem

   SUBROUTINE qes_write_reciprocal_lattice(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(reciprocal_lattice_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%b1)
     CALL xml_addCharacters(xp, obj%b2)
     CALL xml_addCharacters(xp, obj%b3)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_reciprocal_lattice

   SUBROUTINE qes_write_electron_control(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(electron_control_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%diagonalization)
     CALL xml_addCharacters(xp, obj%mixing_mode)
     CALL xml_addCharacters(xp, obj%mixing_beta)
     CALL xml_addCharacters(xp, obj%conv_thr)
     CALL xml_addCharacters(xp, obj%mixing_ndim)
     CALL xml_addCharacters(xp, obj%max_nstep)
     CALL xml_addCharacters(xp, obj%real_space_q)
     CALL xml_addCharacters(xp, obj%tq_smoothing)
     CALL xml_addCharacters(xp, obj%tbeta_smoothing)
     CALL xml_addCharacters(xp, obj%diago_thr_init)
     CALL xml_addCharacters(xp, obj%diago_full_acc)
     IF (obj%diago_cg_maxiter_ispresent) THEN
        CALL xml_addCharacters(xp, obj%diago_cg_maxiter)
     END IF
     IF (obj%diago_david_ndim_ispresent) THEN
        CALL xml_addCharacters(xp, obj%diago_david_ndim)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_electron_control

   SUBROUTINE qes_write_k_points_IBZ(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(k_points_IBZ_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%monkhorst_pack_ispresent) THEN
        CALL qes_write_monkhorst_pack (xp, obj%monkhorst_pack)
     END IF
     IF (obj%nk_ispresent) THEN
        CALL xml_addCharacters(xp, obj%nk)
     END IF
     IF (obj%k_point_ispresent) THEN
        DO i = 1, obj%ndim_k_point
           CALL qes_write_k_point(xp, obj%k_point(i) )
        END DO
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_k_points_IBZ

   SUBROUTINE qes_write_monkhorst_pack(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(monkhorst_pack_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'nk1', obj%nk1 )
     CALL xml_addAttribute(xp, 'nk2', obj%nk2 )
     CALL xml_addAttribute(xp, 'nk3', obj%nk3 )
     CALL xml_addAttribute(xp, 'k1', obj%k1 )
     CALL xml_addAttribute(xp, 'k2', obj%k2 )
     CALL xml_addAttribute(xp, 'k3', obj%k3 )
     CALL xml_AddCharacters(xp, obj%monkhorst_pack)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_monkhorst_pack

   SUBROUTINE qes_write_k_point(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(k_point_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%weight_ispresent) CALL xml_addAttribute(xp, 'weight', obj%weight )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     CALL xml_AddCharacters(xp, obj%k_point)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_k_point

   SUBROUTINE qes_write_ion_control(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(ion_control_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%ion_dynamics)
     IF (obj%upscale_ispresent) THEN
        CALL xml_addCharacters(xp, obj%upscale)
     END IF
     IF (obj%remove_rigid_rot_ispresent) THEN
        CALL xml_addCharacters(xp, obj%remove_rigid_rot)
     END IF
     IF (obj%refold_pos_ispresent) THEN
        CALL xml_addCharacters(xp, obj%refold_pos)
     END IF
     IF (obj%bfgs_ispresent) THEN
        CALL qes_write_bfgs (xp, obj%bfgs)
     END IF
     IF (obj%md_ispresent) THEN
        CALL qes_write_md (xp, obj%md)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_ion_control

   SUBROUTINE qes_write_bfgs(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(bfgs_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%ndim)
     CALL xml_addCharacters(xp, obj%trust_radius_min)
     CALL xml_addCharacters(xp, obj%trust_radius_max)
     CALL xml_addCharacters(xp, obj%trust_radius_init)
     CALL xml_addCharacters(xp, obj%w1)
     CALL xml_addCharacters(xp, obj%w2)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_bfgs

   SUBROUTINE qes_write_md(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(md_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%pot_extrapolation)
     CALL xml_addCharacters(xp, obj%wfc_extrapolation)
     CALL xml_addCharacters(xp, obj%ion_temperature)
     CALL xml_addCharacters(xp, obj%timestep)
     CALL xml_addCharacters(xp, obj%tempw)
     CALL xml_addCharacters(xp, obj%tolp)
     CALL xml_addCharacters(xp, obj%deltaT)
     CALL xml_addCharacters(xp, obj%nraise)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_md

   SUBROUTINE qes_write_cell_control(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cell_control_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%cell_dynamics)
     CALL xml_addCharacters(xp, obj%pressure)
     IF (obj%wmass_ispresent) THEN
        CALL xml_addCharacters(xp, obj%wmass)
     END IF
     IF (obj%cell_factor_ispresent) THEN
        CALL xml_addCharacters(xp, obj%cell_factor)
     END IF
     IF (obj%fix_volume_ispresent) THEN
        CALL xml_addCharacters(xp, obj%fix_volume)
     END IF
     IF (obj%fix_area_ispresent) THEN
        CALL xml_addCharacters(xp, obj%fix_area)
     END IF
     IF (obj%isotropic_ispresent) THEN
        CALL xml_addCharacters(xp, obj%isotropic)
     END IF
     IF (obj%free_cell_ispresent) THEN
        CALL qes_write_integerMatrix (xp, obj%free_cell)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cell_control

   SUBROUTINE qes_write_symmetry_flags(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(symmetry_flags_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%nosym)
     CALL xml_addCharacters(xp, obj%nosym_evc)
     CALL xml_addCharacters(xp, obj%noinv)
     CALL xml_addCharacters(xp, obj%no_t_rev)
     CALL xml_addCharacters(xp, obj%force_symmorphic)
     CALL xml_addCharacters(xp, obj%use_all_frac)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_symmetry_flags

   SUBROUTINE qes_write_boundary_conditions(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(boundary_conditions_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%assume_isolated)
     IF (obj%esm_ispresent) THEN
        CALL qes_write_esm (xp, obj%esm)
     END IF
     IF (obj%fcp_opt_ispresent) THEN
        CALL xml_addCharacters(xp, obj%fcp_opt)
     END IF
     IF (obj%fcp_mu_ispresent) THEN
        CALL xml_addCharacters(xp, obj%fcp_mu)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_boundary_conditions

   SUBROUTINE qes_write_esm(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(esm_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%bc)
     CALL xml_addCharacters(xp, obj%nfit)
     CALL xml_addCharacters(xp, obj%w)
     CALL xml_addCharacters(xp, obj%efield)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_esm

   SUBROUTINE qes_write_ekin_functional(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(ekin_functional_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%ecfixed)
     CALL xml_addCharacters(xp, obj%qcutz)
     CALL xml_addCharacters(xp, obj%q2sigma)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_ekin_functional

   SUBROUTINE qes_write_spin_constraints(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(spin_constraints_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%spin_constraints)
     CALL xml_addCharacters(xp, obj%lagrange_multiplier)
     IF (obj%target_magnetization_ispresent) THEN
        CALL xml_addCharacters(xp, obj%target_magnetization)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_spin_constraints

   SUBROUTINE qes_write_electric_field(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(electric_field_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%electric_potential)
     IF (obj%dipole_correction_ispresent) THEN
        CALL xml_addCharacters(xp, obj%dipole_correction)
     END IF
     IF (obj%gate_settings_ispresent) THEN
        CALL qes_write_gate_settings (xp, obj%gate_settings)
     END IF
     IF (obj%electric_field_direction_ispresent) THEN
        CALL xml_addCharacters(xp, obj%electric_field_direction)
     END IF
     IF (obj%potential_max_position_ispresent) THEN
        CALL xml_addCharacters(xp, obj%potential_max_position)
     END IF
     IF (obj%potential_decrease_width_ispresent) THEN
        CALL xml_addCharacters(xp, obj%potential_decrease_width)
     END IF
     IF (obj%electric_field_amplitude_ispresent) THEN
        CALL xml_addCharacters(xp, obj%electric_field_amplitude)
     END IF
     IF (obj%electric_field_vector_ispresent) THEN
        CALL xml_addCharacters(xp, obj%electric_field_vector)
     END IF
     IF (obj%nk_per_string_ispresent) THEN
        CALL xml_addCharacters(xp, obj%nk_per_string)
     END IF
     IF (obj%n_berry_cycles_ispresent) THEN
        CALL xml_addCharacters(xp, obj%n_berry_cycles)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_electric_field

   SUBROUTINE qes_write_gate_settings(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(gate_settings_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%use_gate)
     IF (obj%zgate_ispresent) THEN
        CALL xml_addCharacters(xp, obj%zgate)
     END IF
     IF (obj%relaxz_ispresent) THEN
        CALL xml_addCharacters(xp, obj%relaxz)
     END IF
     IF (obj%block_ispresent) THEN
        CALL xml_addCharacters(xp, obj%block)
     END IF
     IF (obj%block_1_ispresent) THEN
        CALL xml_addCharacters(xp, obj%block_1)
     END IF
     IF (obj%block_2_ispresent) THEN
        CALL xml_addCharacters(xp, obj%block_2)
     END IF
     IF (obj%block_height_ispresent) THEN
        CALL xml_addCharacters(xp, obj%block_height)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_gate_settings

   SUBROUTINE qes_write_atomic_constraints(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_constraints_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%num_of_constraints)
     CALL xml_addCharacters(xp, obj%tolerance)
     DO i = 1, obj%ndim_atomic_constraint
        CALL qes_write_atomic_constraint(xp, obj%atomic_constraint(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_constraints

   SUBROUTINE qes_write_atomic_constraint(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_constraint_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%constr_parms)
     CALL xml_addCharacters(xp, obj%constr_type)
     CALL xml_addCharacters(xp, obj%constr_target)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_constraint

   SUBROUTINE qes_write_inputOccupations(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(inputOccupations_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'ispin', obj%ispin )
     CALL xml_addAttribute(xp, 'spin_factor', obj%spin_factor )
     CALL xml_addAttribute(xp, 'size', obj%size )
     CALL xml_AddCharacters(xp, obj%inputOccupations, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_inputOccupations

   SUBROUTINE qes_write_outputElectricField(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(outputElectricField_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%BerryPhase_ispresent) THEN
        CALL qes_write_BerryPhaseOutput (xp, obj%BerryPhase)
     END IF
     IF (obj%finiteElectricFieldInfo_ispresent) THEN
        CALL qes_write_finiteFieldOut (xp, obj%finiteElectricFieldInfo)
     END IF
     IF (obj%dipoleInfo_ispresent) THEN
        CALL qes_write_dipoleOutput (xp, obj%dipoleInfo)
     END IF
     IF (obj%gateInfo_ispresent) THEN
        CALL qes_write_gateInfo (xp, obj%gateInfo)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_outputElectricField

   SUBROUTINE qes_write_BerryPhaseOutput(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(BerryPhaseOutput_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_polarization (xp, obj%totalPolarization)
     CALL qes_write_phase (xp, obj%totalPhase)
     DO i = 1, obj%ndim_ionicPolarization
        CALL qes_write_ionicPolarization(xp, obj%ionicPolarization(i) )
     END DO
     DO i = 1, obj%ndim_electronicPolarization
        CALL qes_write_electronicPolarization(xp, obj%electronicPolarization(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_BerryPhaseOutput

   SUBROUTINE qes_write_dipoleOutput(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(dipoleOutput_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%idir)
     CALL qes_write_scalarQuantity (xp, obj%dipole)
     CALL qes_write_scalarQuantity (xp, obj%ion_dipole)
     CALL qes_write_scalarQuantity (xp, obj%elec_dipole)
     CALL qes_write_scalarQuantity (xp, obj%dipoleField)
     CALL qes_write_scalarQuantity (xp, obj%potentialAmp)
     CALL qes_write_scalarQuantity (xp, obj%totalLength)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_dipoleOutput

   SUBROUTINE qes_write_finiteFieldOut(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(finiteFieldOut_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%electronicDipole)
     CALL xml_addCharacters(xp, obj%ionicDipole)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_finiteFieldOut

   SUBROUTINE qes_write_polarization(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(polarization_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_scalarQuantity (xp, obj%polarization)
     CALL xml_addCharacters(xp, obj%modulus)
     CALL xml_addCharacters(xp, obj%direction)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_polarization

   SUBROUTINE qes_write_ionicPolarization(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(ionicPolarization_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_atom (xp, obj%ion)
     CALL xml_addCharacters(xp, obj%charge)
     CALL qes_write_phase (xp, obj%phase)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_ionicPolarization

   SUBROUTINE qes_write_electronicPolarization(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(electronicPolarization_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_k_point (xp, obj%firstKeyPoint)
     IF (obj%spin_ispresent) THEN
        CALL xml_addCharacters(xp, obj%spin)
     END IF
     CALL qes_write_phase (xp, obj%phase)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_electronicPolarization

   SUBROUTINE qes_write_phase(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(phase_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%ionic_ispresent) CALL xml_addAttribute(xp, 'ionic', obj%ionic )
     IF (obj%electronic_ispresent) CALL xml_addAttribute(xp, 'electronic', obj%electronic )
     IF (obj%modulus_ispresent) CALL xml_addAttribute(xp, 'modulus', TRIM(obj%modulus) )
     CALL xml_AddCharacters(xp, obj%phase)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_phase

   SUBROUTINE qes_write_gateInfo(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(gateInfo_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%pot_prefactor)
     CALL xml_addCharacters(xp, obj%gate_zpos)
     CALL xml_addCharacters(xp, obj%gate_gate_term)
     CALL xml_addCharacters(xp, obj%gatefieldEnergy)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_gateInfo

   SUBROUTINE qes_write_convergence_info(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(convergence_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_scf_conv (xp, obj%scf_conv)
     IF (obj%opt_conv_ispresent) THEN
        CALL qes_write_opt_conv (xp, obj%opt_conv)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_convergence_info

   SUBROUTINE qes_write_scf_conv(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(scf_conv_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%n_scf_steps)
     CALL xml_addCharacters(xp, obj%scf_error)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_scf_conv

   SUBROUTINE qes_write_opt_conv(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(opt_conv_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%n_opt_steps)
     CALL xml_addCharacters(xp, obj%grad_norm)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_opt_conv

   SUBROUTINE qes_write_algorithmic_info(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(algorithmic_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%real_space_q)
     CALL xml_addCharacters(xp, obj%uspp)
     CALL xml_addCharacters(xp, obj%paw)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_algorithmic_info

   SUBROUTINE qes_write_symmetries(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(symmetries_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%nsym)
     CALL xml_addCharacters(xp, obj%nrot)
     CALL xml_addCharacters(xp, obj%space_group)
     DO i = 1, obj%ndim_symmetry
        CALL qes_write_symmetry(xp, obj%symmetry(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_symmetries

   SUBROUTINE qes_write_symmetry(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(symmetry_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_info (xp, obj%info)
     CALL qes_write_matrix (xp, obj%rotation)
     IF (obj%fractional_translation_ispresent) THEN
        CALL xml_addCharacters(xp, obj%fractional_translation)
     END IF
     IF (obj%equivalent_atoms_ispresent) THEN
        CALL qes_write_equivalent_atoms (xp, obj%equivalent_atoms)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_symmetry

   SUBROUTINE qes_write_equivalent_atoms(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(equivalent_atoms_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'nat', obj%nat )
     CALL xml_addAttribute(xp, 'size', obj%size )
     CALL xml_AddCharacters(xp, obj%equivalent_atoms)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_equivalent_atoms

   SUBROUTINE qes_write_info(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%name_ispresent) CALL xml_addAttribute(xp, 'name', TRIM(obj%name) )
     IF (obj%class_ispresent) CALL xml_addAttribute(xp, 'class', TRIM(obj%class) )
     IF (obj%time_reversal_ispresent) CALL xml_addAttribute(xp, 'time_reversal', obj%time_reversal )
     CALL xml_AddCharacters(xp, obj%info)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_info

   SUBROUTINE qes_write_outputPBC(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(outputPBC_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%assume_isolated)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_outputPBC

   SUBROUTINE qes_write_magnetization(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(magnetization_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%lsda)
     CALL xml_addCharacters(xp, obj%noncolin)
     CALL xml_addCharacters(xp, obj%spinorbit)
     CALL xml_addCharacters(xp, obj%total)
     CALL xml_addCharacters(xp, obj%absolute)
     CALL xml_addCharacters(xp, obj%do_magnetization)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_magnetization

   SUBROUTINE qes_write_total_energy(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(total_energy_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%etot)
     IF (obj%eband_ispresent) THEN
        CALL xml_addCharacters(xp, obj%eband)
     END IF
     IF (obj%ehart_ispresent) THEN
        CALL xml_addCharacters(xp, obj%ehart)
     END IF
     IF (obj%vtxc_ispresent) THEN
        CALL xml_addCharacters(xp, obj%vtxc)
     END IF
     IF (obj%etxc_ispresent) THEN
        CALL xml_addCharacters(xp, obj%etxc)
     END IF
     IF (obj%ewald_ispresent) THEN
        CALL xml_addCharacters(xp, obj%ewald)
     END IF
     IF (obj%demet_ispresent) THEN
        CALL xml_addCharacters(xp, obj%demet)
     END IF
     IF (obj%efieldcorr_ispresent) THEN
        CALL xml_addCharacters(xp, obj%efieldcorr)
     END IF
     IF (obj%potentiostat_contr_ispresent) THEN
        CALL xml_addCharacters(xp, obj%potentiostat_contr)
     END IF
     IF (obj%gatefield_contr_ispresent) THEN
        CALL xml_addCharacters(xp, obj%gatefield_contr)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_total_energy

   SUBROUTINE qes_write_band_structure(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(band_structure_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addCharacters(xp, obj%lsda)
     CALL xml_addCharacters(xp, obj%noncolin)
     CALL xml_addCharacters(xp, obj%spinorbit)
     CALL xml_addCharacters(xp, obj%nbnd)
     IF (obj%nbnd_up_ispresent) THEN
        CALL xml_addCharacters(xp, obj%nbnd_up)
     END IF
     IF (obj%nbnd_dw_ispresent) THEN
        CALL xml_addCharacters(xp, obj%nbnd_dw)
     END IF
     CALL xml_addCharacters(xp, obj%nelec)
     IF (obj%num_of_atomic_wfc_ispresent) THEN
        CALL xml_addCharacters(xp, obj%num_of_atomic_wfc)
     END IF
     CALL xml_addCharacters(xp, obj%wf_collected)
     IF (obj%fermi_energy_ispresent) THEN
        CALL xml_addCharacters(xp, obj%fermi_energy)
     END IF
     IF (obj%highestOccupiedLevel_ispresent) THEN
        CALL xml_addCharacters(xp, obj%highestOccupiedLevel)
     END IF
     IF (obj%two_fermi_energies_ispresent) THEN
        CALL xml_addCharacters(xp, obj%two_fermi_energies)
     END IF
     CALL qes_write_k_points_IBZ (xp, obj%starting_k_points)
     CALL xml_addCharacters(xp, obj%nks)
     CALL qes_write_occupations (xp, obj%occupations_kind)
     IF (obj%smearing_ispresent) THEN
        CALL qes_write_smearing (xp, obj%smearing)
     END IF
     DO i = 1, obj%ndim_ks_energies
        CALL qes_write_ks_energies(xp, obj%ks_energies(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_band_structure

   SUBROUTINE qes_write_ks_energies(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(ks_energies_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_k_point (xp, obj%k_point)
     CALL xml_addCharacters(xp, obj%npw)
     CALL qes_write_vector (xp, obj%eigenvalues)
     CALL qes_write_vector (xp, obj%occupations)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_ks_energies

   SUBROUTINE qes_write_closed(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(closed_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'DATE', TRIM(obj%DATE) )
     CALL xml_addAttribute(xp, 'TIME', TRIM(obj%TIME) )
     CALL xml_AddCharacters(xp, obj%closed)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_closed

   SUBROUTINE qes_write_vector(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(vector_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'size', obj%size )
     CALL xml_AddCharacters(xp, obj%vector, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_vector

   SUBROUTINE qes_write_integerVector(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(integerVector_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'size', obj%size )
     CALL xml_AddCharacters(xp, obj%integerVector)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_integerVector

   SUBROUTINE qes_write_matrix(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(matrix_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'rank', obj%rank )
     CALL xml_addAttribute(xp, 'dims', obj%dims )
     CALL xml_addAttribute(xp, 'order', TRIM(obj%order) )
       CALL xml_addNewLine(xp)
        DO i = 1, obj%dims(2)
           CALL xml_AddCharacters(xp, obj%matrix((i-1)*obj%dims(1)+1: i*obj%dims(1)), fmt ='s16')
           CALL xml_addNewLine(xp)
        END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_matrix

   SUBROUTINE qes_write_integerMatrix(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(integerMatrix_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'rank', obj%rank )
     CALL xml_addAttribute(xp, 'dims', obj%dims )
     CALL xml_addAttribute(xp, 'order', TRIM(obj%order) )
        CALL xml_addNewLine(xp)
        DO i = 1, obj%dims(2)
           CALL xml_AddCharacters(xp, obj%integerMatrix((i-1)*obj%dims(1)+1: i*obj%dims(1)) )
           CALL xml_addNewLine(xp)
        END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_integerMatrix

   SUBROUTINE qes_write_scalarQuantity(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(scalarQuantity_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'Units', TRIM(obj%Units) )
     CALL xml_AddCharacters(xp, obj%scalarQuantity)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_scalarQuantity

  !
END MODULE qes_write_module