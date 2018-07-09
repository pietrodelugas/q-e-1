!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qes_libs_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes-1.0
  !
  USE qes_init_module
  USE qes_reset_module
  USE qes_read_module
  USE qes_write_module
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: qes_init, qes_bcast, qes_reset, qes_read, qes_write
  !
  !INTERFACE qes_init
  !
  !  MODULE PROCEDURE qes_init_general_info
  !
  !  MODULE PROCEDURE qes_init_parallel_info
  !
  !  MODULE PROCEDURE qes_init_input
  !
  !  MODULE PROCEDURE qes_init_step
  !
  !  MODULE PROCEDURE qes_init_output
  !
  !  MODULE PROCEDURE qes_init_control_variables
  !
  !  MODULE PROCEDURE qes_init_xml_format
  !
  !  MODULE PROCEDURE qes_init_creator
  !
  !  MODULE PROCEDURE qes_init_created
  !
  !  MODULE PROCEDURE qes_init_atomic_species
  !
  !  MODULE PROCEDURE qes_init_species
  !
  !  MODULE PROCEDURE qes_init_atomic_structure
  !
  !  MODULE PROCEDURE qes_init_atomic_positions
  !
  !  MODULE PROCEDURE qes_init_atom
  !
  !  MODULE PROCEDURE qes_init_wyckoff_positions
  !
  !  MODULE PROCEDURE qes_init_cell
  !
  !  MODULE PROCEDURE qes_init_dft
  !
  !  MODULE PROCEDURE qes_init_hybrid
  !
  !  MODULE PROCEDURE qes_init_qpoint_grid
  !
  !  MODULE PROCEDURE qes_init_dftU
  !
  !  MODULE PROCEDURE qes_init_HubbardCommon
  !
  !  MODULE PROCEDURE qes_init_HubbardJ
  !
  !  MODULE PROCEDURE qes_init_starting_ns
  !
  !  MODULE PROCEDURE qes_init_Hubbard_ns
  !
  !  MODULE PROCEDURE qes_init_vdW
  !
  !  MODULE PROCEDURE qes_init_spin
  !
  !  MODULE PROCEDURE qes_init_bands
  !
  !  MODULE PROCEDURE qes_init_smearing
  !
  !  MODULE PROCEDURE qes_init_occupations
  !
  !  MODULE PROCEDURE qes_init_basis
  !
  !  MODULE PROCEDURE qes_init_basis_set
  !
  !  MODULE PROCEDURE qes_init_basisSetItem
  !
  !  MODULE PROCEDURE qes_init_reciprocal_lattice
  !
  !  MODULE PROCEDURE qes_init_electron_control
  !
  !  MODULE PROCEDURE qes_init_k_points_IBZ
  !
  !  MODULE PROCEDURE qes_init_monkhorst_pack
  !
  !  MODULE PROCEDURE qes_init_k_point
  !
  !  MODULE PROCEDURE qes_init_ion_control
  !
  !  MODULE PROCEDURE qes_init_bfgs
  !
  !  MODULE PROCEDURE qes_init_md
  !
  !  MODULE PROCEDURE qes_init_cell_control
  !
  !  MODULE PROCEDURE qes_init_symmetry_flags
  !
  !  MODULE PROCEDURE qes_init_boundary_conditions
  !
  !  MODULE PROCEDURE qes_init_esm
  !
  !  MODULE PROCEDURE qes_init_ekin_functional
  !
  !  MODULE PROCEDURE qes_init_spin_constraints
  !
  !  MODULE PROCEDURE qes_init_electric_field
  !
  !  MODULE PROCEDURE qes_init_gate_settings
  !
  !  MODULE PROCEDURE qes_init_atomic_constraints
  !
  !  MODULE PROCEDURE qes_init_atomic_constraint
  !
  !  MODULE PROCEDURE qes_init_inputOccupations
  !
  !  MODULE PROCEDURE qes_init_outputElectricField
  !
  !  MODULE PROCEDURE qes_init_BerryPhaseOutput
  !
  !  MODULE PROCEDURE qes_init_dipoleOutput
  !
  !  MODULE PROCEDURE qes_init_finiteFieldOut
  !
  !  MODULE PROCEDURE qes_init_polarization
  !
  !  MODULE PROCEDURE qes_init_ionicPolarization
  !
  !  MODULE PROCEDURE qes_init_electronicPolarization
  !
  !  MODULE PROCEDURE qes_init_phase
  !
  !  MODULE PROCEDURE qes_init_gateInfo
  !
  !  MODULE PROCEDURE qes_init_convergence_info
  !
  !  MODULE PROCEDURE qes_init_scf_conv
  !
  !  MODULE PROCEDURE qes_init_opt_conv
  !
  !  MODULE PROCEDURE qes_init_algorithmic_info
  !
  !  MODULE PROCEDURE qes_init_symmetries
  !
  !  MODULE PROCEDURE qes_init_symmetry
  !
  !  MODULE PROCEDURE qes_init_equivalent_atoms
  !
  !  MODULE PROCEDURE qes_init_info
  !
  !  MODULE PROCEDURE qes_init_outputPBC
  !
  !  MODULE PROCEDURE qes_init_magnetization
  !
  !  MODULE PROCEDURE qes_init_total_energy
  !
  !  MODULE PROCEDURE qes_init_band_structure
  !
  !  MODULE PROCEDURE qes_init_ks_energies
  !
  !  MODULE PROCEDURE qes_init_closed
  !
  !  MODULE PROCEDURE qes_init_vector
  !
  !  MODULE PROCEDURE qes_init_integerVector
  !
  !  MODULE PROCEDURE qes_init_matrix_1, qes_init_matrix_2, qes_init_matrix_3
  !
  !  MODULE PROCEDURE qes_init_integerMatrix_1, qes_init_integerMatrix_2, qes_init_integerMatrix_3
  !
  !  MODULE PROCEDURE qes_init_scalarQuantity
  !
  !END INTERFACE qes_init
  !
  !INTERFACE qes_bcast
  !
  !  MODULE PROCEDURE qes_bcast_general_info
  !
  !  MODULE PROCEDURE qes_bcast_parallel_info
  !
  !  MODULE PROCEDURE qes_bcast_input
  !
  !  MODULE PROCEDURE qes_bcast_step
  !
  !  MODULE PROCEDURE qes_bcast_output
  !
  !  MODULE PROCEDURE qes_bcast_control_variables
  !
  !  MODULE PROCEDURE qes_bcast_xml_format
  !
  !  MODULE PROCEDURE qes_bcast_creator
  !
  !  MODULE PROCEDURE qes_bcast_created
  !
  !  MODULE PROCEDURE qes_bcast_atomic_species
  !
  !  MODULE PROCEDURE qes_bcast_species
  !
  !  MODULE PROCEDURE qes_bcast_atomic_structure
  !
  !  MODULE PROCEDURE qes_bcast_atomic_positions
  !
  !  MODULE PROCEDURE qes_bcast_atom
  !
  !  MODULE PROCEDURE qes_bcast_wyckoff_positions
  !
  !  MODULE PROCEDURE qes_bcast_cell
  !
  !  MODULE PROCEDURE qes_bcast_dft
  !
  !  MODULE PROCEDURE qes_bcast_hybrid
  !
  !  MODULE PROCEDURE qes_bcast_qpoint_grid
  !
  !  MODULE PROCEDURE qes_bcast_dftU
  !
  !  MODULE PROCEDURE qes_bcast_HubbardCommon
  !
  !  MODULE PROCEDURE qes_bcast_HubbardJ
  !
  !  MODULE PROCEDURE qes_bcast_starting_ns
  !
  !  MODULE PROCEDURE qes_bcast_Hubbard_ns
  !
  !  MODULE PROCEDURE qes_bcast_vdW
  !
  !  MODULE PROCEDURE qes_bcast_spin
  !
  !  MODULE PROCEDURE qes_bcast_bands
  !
  !  MODULE PROCEDURE qes_bcast_smearing
  !
  !  MODULE PROCEDURE qes_bcast_occupations
  !
  !  MODULE PROCEDURE qes_bcast_basis
  !
  !  MODULE PROCEDURE qes_bcast_basis_set
  !
  !  MODULE PROCEDURE qes_bcast_basisSetItem
  !
  !  MODULE PROCEDURE qes_bcast_reciprocal_lattice
  !
  !  MODULE PROCEDURE qes_bcast_electron_control
  !
  !  MODULE PROCEDURE qes_bcast_k_points_IBZ
  !
  !  MODULE PROCEDURE qes_bcast_monkhorst_pack
  !
  !  MODULE PROCEDURE qes_bcast_k_point
  !
  !  MODULE PROCEDURE qes_bcast_ion_control
  !
  !  MODULE PROCEDURE qes_bcast_bfgs
  !
  !  MODULE PROCEDURE qes_bcast_md
  !
  !  MODULE PROCEDURE qes_bcast_cell_control
  !
  !  MODULE PROCEDURE qes_bcast_symmetry_flags
  !
  !  MODULE PROCEDURE qes_bcast_boundary_conditions
  !
  !  MODULE PROCEDURE qes_bcast_esm
  !
  !  MODULE PROCEDURE qes_bcast_ekin_functional
  !
  !  MODULE PROCEDURE qes_bcast_spin_constraints
  !
  !  MODULE PROCEDURE qes_bcast_electric_field
  !
  !  MODULE PROCEDURE qes_bcast_gate_settings
  !
  !  MODULE PROCEDURE qes_bcast_atomic_constraints
  !
  !  MODULE PROCEDURE qes_bcast_atomic_constraint
  !
  !  MODULE PROCEDURE qes_bcast_inputOccupations
  !
  !  MODULE PROCEDURE qes_bcast_outputElectricField
  !
  !  MODULE PROCEDURE qes_bcast_BerryPhaseOutput
  !
  !  MODULE PROCEDURE qes_bcast_dipoleOutput
  !
  !  MODULE PROCEDURE qes_bcast_finiteFieldOut
  !
  !  MODULE PROCEDURE qes_bcast_polarization
  !
  !  MODULE PROCEDURE qes_bcast_ionicPolarization
  !
  !  MODULE PROCEDURE qes_bcast_electronicPolarization
  !
  !  MODULE PROCEDURE qes_bcast_phase
  !
  !  MODULE PROCEDURE qes_bcast_gateInfo
  !
  !  MODULE PROCEDURE qes_bcast_convergence_info
  !
  !  MODULE PROCEDURE qes_bcast_scf_conv
  !
  !  MODULE PROCEDURE qes_bcast_opt_conv
  !
  !  MODULE PROCEDURE qes_bcast_algorithmic_info
  !
  !  MODULE PROCEDURE qes_bcast_symmetries
  !
  !  MODULE PROCEDURE qes_bcast_symmetry
  !
  !  MODULE PROCEDURE qes_bcast_equivalent_atoms
  !
  !  MODULE PROCEDURE qes_bcast_info
  !
  !  MODULE PROCEDURE qes_bcast_outputPBC
  !
  !  MODULE PROCEDURE qes_bcast_magnetization
  !
  !  MODULE PROCEDURE qes_bcast_total_energy
  !
  !  MODULE PROCEDURE qes_bcast_band_structure
  !
  !  MODULE PROCEDURE qes_bcast_ks_energies
  !
  !  MODULE PROCEDURE qes_bcast_closed
  !
  !  MODULE PROCEDURE qes_bcast_vector
  !
  !  MODULE PROCEDURE qes_bcast_integerVector
  !
  !  MODULE PROCEDURE qes_bcast_matrix
  !
  !  MODULE PROCEDURE qes_bcast_integerMatrix
  !
  !  MODULE PROCEDURE qes_bcast_scalarQuantity
  !
  !END INTERFACE qes_bcast
  !
  !INTERFACE qes_read
  !
  !  MODULE PROCEDURE qes_read_general_info
  !
  !  MODULE PROCEDURE qes_read_parallel_info
  !
  !  MODULE PROCEDURE qes_read_input
  !
  !  MODULE PROCEDURE qes_read_step
  !
  !  MODULE PROCEDURE qes_read_output
  !
  !  MODULE PROCEDURE qes_read_control_variables
  !
  !  MODULE PROCEDURE qes_read_xml_format
  !
  !  MODULE PROCEDURE qes_read_creator
  !
  !  MODULE PROCEDURE qes_read_created
  !
  !  MODULE PROCEDURE qes_read_atomic_species
  !
  !  MODULE PROCEDURE qes_read_species
  !
  !  MODULE PROCEDURE qes_read_atomic_structure
  !
  !  MODULE PROCEDURE qes_read_atomic_positions
  !
  !  MODULE PROCEDURE qes_read_atom
  !
  !  MODULE PROCEDURE qes_read_wyckoff_positions
  !
  !  MODULE PROCEDURE qes_read_cell
  !
  !  MODULE PROCEDURE qes_read_dft
  !
  !  MODULE PROCEDURE qes_read_hybrid
  !
  !  MODULE PROCEDURE qes_read_qpoint_grid
  !
  !  MODULE PROCEDURE qes_read_dftU
  !
  !  MODULE PROCEDURE qes_read_HubbardCommon
  !
  !  MODULE PROCEDURE qes_read_HubbardJ
  !
  !  MODULE PROCEDURE qes_read_starting_ns
  !
  !  MODULE PROCEDURE qes_read_Hubbard_ns
  !
  !  MODULE PROCEDURE qes_read_vdW
  !
  !  MODULE PROCEDURE qes_read_spin
  !
  !  MODULE PROCEDURE qes_read_bands
  !
  !  MODULE PROCEDURE qes_read_smearing
  !
  !  MODULE PROCEDURE qes_read_occupations
  !
  !  MODULE PROCEDURE qes_read_basis
  !
  !  MODULE PROCEDURE qes_read_basis_set
  !
  !  MODULE PROCEDURE qes_read_basisSetItem
  !
  !  MODULE PROCEDURE qes_read_reciprocal_lattice
  !
  !  MODULE PROCEDURE qes_read_electron_control
  !
  !  MODULE PROCEDURE qes_read_k_points_IBZ
  !
  !  MODULE PROCEDURE qes_read_monkhorst_pack
  !
  !  MODULE PROCEDURE qes_read_k_point
  !
  !  MODULE PROCEDURE qes_read_ion_control
  !
  !  MODULE PROCEDURE qes_read_bfgs
  !
  !  MODULE PROCEDURE qes_read_md
  !
  !  MODULE PROCEDURE qes_read_cell_control
  !
  !  MODULE PROCEDURE qes_read_symmetry_flags
  !
  !  MODULE PROCEDURE qes_read_boundary_conditions
  !
  !  MODULE PROCEDURE qes_read_esm
  !
  !  MODULE PROCEDURE qes_read_ekin_functional
  !
  !  MODULE PROCEDURE qes_read_spin_constraints
  !
  !  MODULE PROCEDURE qes_read_electric_field
  !
  !  MODULE PROCEDURE qes_read_gate_settings
  !
  !  MODULE PROCEDURE qes_read_atomic_constraints
  !
  !  MODULE PROCEDURE qes_read_atomic_constraint
  !
  !  MODULE PROCEDURE qes_read_inputOccupations
  !
  !  MODULE PROCEDURE qes_read_outputElectricField
  !
  !  MODULE PROCEDURE qes_read_BerryPhaseOutput
  !
  !  MODULE PROCEDURE qes_read_dipoleOutput
  !
  !  MODULE PROCEDURE qes_read_finiteFieldOut
  !
  !  MODULE PROCEDURE qes_read_polarization
  !
  !  MODULE PROCEDURE qes_read_ionicPolarization
  !
  !  MODULE PROCEDURE qes_read_electronicPolarization
  !
  !  MODULE PROCEDURE qes_read_phase
  !
  !  MODULE PROCEDURE qes_read_gateInfo
  !
  !  MODULE PROCEDURE qes_read_convergence_info
  !
  !  MODULE PROCEDURE qes_read_scf_conv
  !
  !  MODULE PROCEDURE qes_read_opt_conv
  !
  !  MODULE PROCEDURE qes_read_algorithmic_info
  !
  !  MODULE PROCEDURE qes_read_symmetries
  !
  !  MODULE PROCEDURE qes_read_symmetry
  !
  !  MODULE PROCEDURE qes_read_equivalent_atoms
  !
  !  MODULE PROCEDURE qes_read_info
  !
  !  MODULE PROCEDURE qes_read_outputPBC
  !
  !  MODULE PROCEDURE qes_read_magnetization
  !
  !  MODULE PROCEDURE qes_read_total_energy
  !
  !  MODULE PROCEDURE qes_read_band_structure
  !
  !  MODULE PROCEDURE qes_read_ks_energies
  !
  !  MODULE PROCEDURE qes_read_closed
  !
  !  MODULE PROCEDURE qes_read_vector
  !
  !  MODULE PROCEDURE qes_read_integerVector
  !
  !  MODULE PROCEDURE qes_read_matrix
  !
  !  MODULE PROCEDURE qes_read_integerMatrix
  !
  !  MODULE PROCEDURE qes_read_scalarQuantity
  !
  !END INTERFACE qes_read
  !
  !INTERFACE qes_reset
  !
  !  MODULE PROCEDURE qes_reset_general_info
  !
  !  MODULE PROCEDURE qes_reset_parallel_info
  !
  !  MODULE PROCEDURE qes_reset_input
  !
  !  MODULE PROCEDURE qes_reset_step
  !
  !  MODULE PROCEDURE qes_reset_output
  !
  !  MODULE PROCEDURE qes_reset_control_variables
  !
  !  MODULE PROCEDURE qes_reset_xml_format
  !
  !  MODULE PROCEDURE qes_reset_creator
  !
  !  MODULE PROCEDURE qes_reset_created
  !
  !  MODULE PROCEDURE qes_reset_atomic_species
  !
  !  MODULE PROCEDURE qes_reset_species
  !
  !  MODULE PROCEDURE qes_reset_atomic_structure
  !
  !  MODULE PROCEDURE qes_reset_atomic_positions
  !
  !  MODULE PROCEDURE qes_reset_atom
  !
  !  MODULE PROCEDURE qes_reset_wyckoff_positions
  !
  !  MODULE PROCEDURE qes_reset_cell
  !
  !  MODULE PROCEDURE qes_reset_dft
  !
  !  MODULE PROCEDURE qes_reset_hybrid
  !
  !  MODULE PROCEDURE qes_reset_qpoint_grid
  !
  !  MODULE PROCEDURE qes_reset_dftU
  !
  !  MODULE PROCEDURE qes_reset_HubbardCommon
  !
  !  MODULE PROCEDURE qes_reset_HubbardJ
  !
  !  MODULE PROCEDURE qes_reset_starting_ns
  !
  !  MODULE PROCEDURE qes_reset_Hubbard_ns
  !
  !  MODULE PROCEDURE qes_reset_vdW
  !
  !  MODULE PROCEDURE qes_reset_spin
  !
  !  MODULE PROCEDURE qes_reset_bands
  !
  !  MODULE PROCEDURE qes_reset_smearing
  !
  !  MODULE PROCEDURE qes_reset_occupations
  !
  !  MODULE PROCEDURE qes_reset_basis
  !
  !  MODULE PROCEDURE qes_reset_basis_set
  !
  !  MODULE PROCEDURE qes_reset_basisSetItem
  !
  !  MODULE PROCEDURE qes_reset_reciprocal_lattice
  !
  !  MODULE PROCEDURE qes_reset_electron_control
  !
  !  MODULE PROCEDURE qes_reset_k_points_IBZ
  !
  !  MODULE PROCEDURE qes_reset_monkhorst_pack
  !
  !  MODULE PROCEDURE qes_reset_k_point
  !
  !  MODULE PROCEDURE qes_reset_ion_control
  !
  !  MODULE PROCEDURE qes_reset_bfgs
  !
  !  MODULE PROCEDURE qes_reset_md
  !
  !  MODULE PROCEDURE qes_reset_cell_control
  !
  !  MODULE PROCEDURE qes_reset_symmetry_flags
  !
  !  MODULE PROCEDURE qes_reset_boundary_conditions
  !
  !  MODULE PROCEDURE qes_reset_esm
  !
  !  MODULE PROCEDURE qes_reset_ekin_functional
  !
  !  MODULE PROCEDURE qes_reset_spin_constraints
  !
  !  MODULE PROCEDURE qes_reset_electric_field
  !
  !  MODULE PROCEDURE qes_reset_gate_settings
  !
  !  MODULE PROCEDURE qes_reset_atomic_constraints
  !
  !  MODULE PROCEDURE qes_reset_atomic_constraint
  !
  !  MODULE PROCEDURE qes_reset_inputOccupations
  !
  !  MODULE PROCEDURE qes_reset_outputElectricField
  !
  !  MODULE PROCEDURE qes_reset_BerryPhaseOutput
  !
  !  MODULE PROCEDURE qes_reset_dipoleOutput
  !
  !  MODULE PROCEDURE qes_reset_finiteFieldOut
  !
  !  MODULE PROCEDURE qes_reset_polarization
  !
  !  MODULE PROCEDURE qes_reset_ionicPolarization
  !
  !  MODULE PROCEDURE qes_reset_electronicPolarization
  !
  !  MODULE PROCEDURE qes_reset_phase
  !
  !  MODULE PROCEDURE qes_reset_gateInfo
  !
  !  MODULE PROCEDURE qes_reset_convergence_info
  !
  !  MODULE PROCEDURE qes_reset_scf_conv
  !
  !  MODULE PROCEDURE qes_reset_opt_conv
  !
  !  MODULE PROCEDURE qes_reset_algorithmic_info
  !
  !  MODULE PROCEDURE qes_reset_symmetries
  !
  !  MODULE PROCEDURE qes_reset_symmetry
  !
  !  MODULE PROCEDURE qes_reset_equivalent_atoms
  !
  !  MODULE PROCEDURE qes_reset_info
  !
  !  MODULE PROCEDURE qes_reset_outputPBC
  !
  !  MODULE PROCEDURE qes_reset_magnetization
  !
  !  MODULE PROCEDURE qes_reset_total_energy
  !
  !  MODULE PROCEDURE qes_reset_band_structure
  !
  !  MODULE PROCEDURE qes_reset_ks_energies
  !
  !  MODULE PROCEDURE qes_reset_closed
  !
  !  MODULE PROCEDURE qes_reset_vector
  !
  !  MODULE PROCEDURE qes_reset_integerVector
  !
  !  MODULE PROCEDURE qes_reset_matrix
  !
  !  MODULE PROCEDURE qes_reset_integerMatrix
  !
  !  MODULE PROCEDURE qes_reset_scalarQuantity
  !
  !END INTERFACE qes_reset
  !
  !INTERFACE qes_write
  !
  !  MODULE PROCEDURE qes_write_general_info
  !
  !  MODULE PROCEDURE qes_write_parallel_info
  !
  !  MODULE PROCEDURE qes_write_input
  !
  !  MODULE PROCEDURE qes_write_step
  !
  !  MODULE PROCEDURE qes_write_output
  !
  !  MODULE PROCEDURE qes_write_control_variables
  !
  !  MODULE PROCEDURE qes_write_xml_format
  !
  !  MODULE PROCEDURE qes_write_creator
  !
  !  MODULE PROCEDURE qes_write_created
  !
  !  MODULE PROCEDURE qes_write_atomic_species
  !
  !  MODULE PROCEDURE qes_write_species
  !
  !  MODULE PROCEDURE qes_write_atomic_structure
  !
  !  MODULE PROCEDURE qes_write_atomic_positions
  !
  !  MODULE PROCEDURE qes_write_atom
  !
  !  MODULE PROCEDURE qes_write_wyckoff_positions
  !
  !  MODULE PROCEDURE qes_write_cell
  !
  !  MODULE PROCEDURE qes_write_dft
  !
  !  MODULE PROCEDURE qes_write_hybrid
  !
  !  MODULE PROCEDURE qes_write_qpoint_grid
  !
  !  MODULE PROCEDURE qes_write_dftU
  !
  !  MODULE PROCEDURE qes_write_HubbardCommon
  !
  !  MODULE PROCEDURE qes_write_HubbardJ
  !
  !  MODULE PROCEDURE qes_write_starting_ns
  !
  !  MODULE PROCEDURE qes_write_Hubbard_ns
  !
  !  MODULE PROCEDURE qes_write_vdW
  !
  !  MODULE PROCEDURE qes_write_spin
  !
  !  MODULE PROCEDURE qes_write_bands
  !
  !  MODULE PROCEDURE qes_write_smearing
  !
  !  MODULE PROCEDURE qes_write_occupations
  !
  !  MODULE PROCEDURE qes_write_basis
  !
  !  MODULE PROCEDURE qes_write_basis_set
  !
  !  MODULE PROCEDURE qes_write_basisSetItem
  !
  !  MODULE PROCEDURE qes_write_reciprocal_lattice
  !
  !  MODULE PROCEDURE qes_write_electron_control
  !
  !  MODULE PROCEDURE qes_write_k_points_IBZ
  !
  !  MODULE PROCEDURE qes_write_monkhorst_pack
  !
  !  MODULE PROCEDURE qes_write_k_point
  !
  !  MODULE PROCEDURE qes_write_ion_control
  !
  !  MODULE PROCEDURE qes_write_bfgs
  !
  !  MODULE PROCEDURE qes_write_md
  !
  !  MODULE PROCEDURE qes_write_cell_control
  !
  !  MODULE PROCEDURE qes_write_symmetry_flags
  !
  !  MODULE PROCEDURE qes_write_boundary_conditions
  !
  !  MODULE PROCEDURE qes_write_esm
  !
  !  MODULE PROCEDURE qes_write_ekin_functional
  !
  !  MODULE PROCEDURE qes_write_spin_constraints
  !
  !  MODULE PROCEDURE qes_write_electric_field
  !
  !  MODULE PROCEDURE qes_write_gate_settings
  !
  !  MODULE PROCEDURE qes_write_atomic_constraints
  !
  !  MODULE PROCEDURE qes_write_atomic_constraint
  !
  !  MODULE PROCEDURE qes_write_inputOccupations
  !
  !  MODULE PROCEDURE qes_write_outputElectricField
  !
  !  MODULE PROCEDURE qes_write_BerryPhaseOutput
  !
  !  MODULE PROCEDURE qes_write_dipoleOutput
  !
  !  MODULE PROCEDURE qes_write_finiteFieldOut
  !
  !  MODULE PROCEDURE qes_write_polarization
  !
  !  MODULE PROCEDURE qes_write_ionicPolarization
  !
  !  MODULE PROCEDURE qes_write_electronicPolarization
  !
  !  MODULE PROCEDURE qes_write_phase
  !
  !  MODULE PROCEDURE qes_write_gateInfo
  !
  !  MODULE PROCEDURE qes_write_convergence_info
  !
  !  MODULE PROCEDURE qes_write_scf_conv
  !
  !  MODULE PROCEDURE qes_write_opt_conv
  !
  !  MODULE PROCEDURE qes_write_algorithmic_info
  !
  !  MODULE PROCEDURE qes_write_symmetries
  !
  !  MODULE PROCEDURE qes_write_symmetry
  !
  !  MODULE PROCEDURE qes_write_equivalent_atoms
  !
  !  MODULE PROCEDURE qes_write_info
  !
  !  MODULE PROCEDURE qes_write_outputPBC
  !
  !  MODULE PROCEDURE qes_write_magnetization
  !
  !  MODULE PROCEDURE qes_write_total_energy
  !
  !  MODULE PROCEDURE qes_write_band_structure
  !
  !  MODULE PROCEDURE qes_write_ks_energies
  !
  !  MODULE PROCEDURE qes_write_closed
  !
  !  MODULE PROCEDURE qes_write_vector
  !
  !  MODULE PROCEDURE qes_write_integerVector
  !
  !  MODULE PROCEDURE qes_write_matrix
  !
  !  MODULE PROCEDURE qes_write_integerMatrix
  !
  !  MODULE PROCEDURE qes_write_scalarQuantity
  !
  !END INTERFACE qes_write
  !
END MODULE qes_libs_module