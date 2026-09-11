module dtm_elem
implicit none

logical :: dtm_pressure_mass_matrix,dtm_only_pressure_mass_matrix
logical :: num_dofs_is_constant ! is true for temperature; false for P2P1
logical :: print_matrix, do_solve, do_dtm_elem

namelist / dtm_elem_nml / dtm_pressure_mass_matrix, dtm_only_pressure_mass_matrix, print_matrix, do_solve

end module dtm_elem

