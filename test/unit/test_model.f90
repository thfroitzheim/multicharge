! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
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

module test_model
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, test_failed
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge_charge, only: get_charges, get_eeq_charges, get_eeqbc_charges
   use multicharge_model, only: mchrg_model_type
   use multicharge_model_cache, only: cache_container
   use multicharge_model_eeqbc, only: eeqbc_model
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   implicit none
   private

   public :: collect_model

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5 * epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

!> Collect all exported unit tests
subroutine collect_model(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("eeq-dadr-mb01", test_eeq_dadr_mb01), &
      & new_unittest("eeq-dadL-mb01", test_eeq_dadL_mb01), &
      & new_unittest("eeq-dbdr-mb01", test_eeq_dbdr_mb01), &
      & new_unittest("eeq-dbdL-mb01", test_eeq_dbdL_mb01), &
      & new_unittest("eeq-charges-mb01", test_eeq_q_mb01), &
      & new_unittest("eeq-charges-mb02", test_eeq_q_mb02), &
      & new_unittest("eeq-charges-efield-mb03", test_eeq_q_efield_mb03), &
      & new_unittest("eeq-charges-actinides", test_eeq_q_actinides), &
      & new_unittest("eeq-energy-mb03", test_eeq_e_mb03), &
      & new_unittest("eeq-energy-mb04", test_eeq_e_mb04), &
      & new_unittest("eeq-gradient-mb05", test_eeq_g_mb05), &
      & new_unittest("eeq-gradient-mb06", test_eeq_g_mb06), &
      & new_unittest("eeq-sigma-mb07", test_eeq_s_mb07), &
      & new_unittest("eeq-sigma-mb08", test_eeq_s_mb08), &
      & new_unittest("eeq-dqdr-mb09", test_eeq_dqdr_mb09), &
      & new_unittest("eeq-dqdr-mb10", test_eeq_dqdr_mb10), &
      & new_unittest("eeq-dqdL-mb11", test_eeq_dqdL_mb11), &
      & new_unittest("eeq-dqdL-mb12", test_eeq_dqdL_mb12), &
      & new_unittest("gradient-h2plus", test_g_h2plus), &
      & new_unittest("eeq-dadr-znooh", test_eeq_dadr_znooh), &
      & new_unittest("eeq-dbdr-znooh", test_eeq_dbdr_znooh), &
      & new_unittest("gradient-znooh", test_g_znooh), &
      & new_unittest("dqdr-znooh", test_dqdr_znooh), &
      & new_unittest("eeqbc-dadr-mb01", test_eeqbc_dadr_mb01), &
      & new_unittest("eeqbc-dadL-mb01", test_eeqbc_dadL_mb01), &
      & new_unittest("eeqbc-dbdr-mb01", test_eeqbc_dbdr_mb01), &
      & new_unittest("eeqbc-dbdL-mb01", test_eeqbc_dbdL_mb01), &
      & new_unittest("eeqbc-dadr-mb05", test_eeqbc_dadr_mb05), &
      & new_unittest("eeqbc-dadL-mb05", test_eeqbc_dadL_mb05), &
      & new_unittest("eeqbc-dbdr-mb05", test_eeqbc_dbdr_mb05), &
      & new_unittest("eeqbc-charges-mb01", test_eeqbc_q_mb01), &
      & new_unittest("eeqbc-charges-mb02", test_eeqbc_q_mb02), &
      & new_unittest("eeqbc-charges-efield-mb03", test_eeqbc_q_efield_mb03), &
      & new_unittest("eeqbc-charges-actinides", test_eeqbc_q_actinides), &
      & new_unittest("eeqbc-energy-mb03", test_eeqbc_e_mb03), &
      & new_unittest("eeqbc-energy-mb04", test_eeqbc_e_mb04), &
      & new_unittest("eeqbc-gradient-mb05", test_eeqbc_g_mb05), &
      & new_unittest("eeqbc-gradient-mb06", test_eeqbc_g_mb06), &
      & new_unittest("eeqbc-sigma-mb07", test_eeqbc_s_mb07), &
      & new_unittest("eeqbc-sigma-mb08", test_eeqbc_s_mb08), &
      & new_unittest("eeqbc-dqdr-mb09", test_eeqbc_dqdr_mb09), &
      & new_unittest("eeqbc-dqdr-mb10", test_eeqbc_dqdr_mb10), &
      & new_unittest("eeqbc-dqdL-mb11", test_eeqbc_dqdL_mb11), &
      & new_unittest("eeqbc-dqdL-mb12", test_eeqbc_dqdL_mb12) &
      & ]

end subroutine collect_model

subroutine test_dadr(error, mol, model, thr_in)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic, jat, kat
   real(wp) :: thr_
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 5.0e-6_wp
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: qloc(:)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: qvec(:), numgrad(:, :, :), amatr1(:, :), amatr2(:, :), amatl1(:, :), amatl2(:, :), numtrace(:, :)
   type(cache_container), allocatable :: cache
   allocate (cache)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate (cn(mol%nat), qloc(mol%nat), amatr1(mol%nat + 1, mol%nat + 1), amatl1(mol%nat + 1, mol%nat + 1), &
      & amatr2(mol%nat + 1, mol%nat + 1), amatl2(mol%nat + 1, mol%nat + 1), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), dqlocdr(3, mol%nat, mol%nat), &
      & dqlocdL(3, 3, mol%nat), dadr(3, mol%nat, mol%nat + 1), dadL(3, 3, mol%nat + 1), &
      & atrace(3, mol%nat), numtrace(3, mol%nat), numgrad(3, mol%nat, mol%nat + 1), qvec(mol%nat))

   ! Obtain the vector of charges
   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%solve(mol, error, cn, qloc, qvec=qvec)
   if (allocated(error)) return

   numgrad = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         amatr1(:, :) = 0.0_wp
         amatr2(:, :) = 0.0_wp
         amatl1(:, :) = 0.0_wp
         amatl2(:, :) = 0.0_wp

         ! First right-hand side (x+h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatr1)

         ! Second right-hand side (x+2h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatr2)

         ! Return to original position before calculating left sides
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step

         ! First left-hand side (x-h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl1)

         ! Second left-hand side (x-2h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl2)

         ! Return to original position
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + 2*step

         do kat = 1, mol%nat
            do jat = 1, mol%nat
               ! Numerical gradient using 4-step central difference formula
               ! f'(x) ≈ [f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)]/(12h)
               numgrad(ic, iat, kat) = numgrad(ic, iat, kat) + &
                  & qvec(jat)*(amatl2(kat, jat) - 8.0_wp*amatl1(kat, jat) + &
                  & 8.0_wp*amatr1(kat, jat) - amatr2(kat, jat))/(12.0_wp*step)
            end do
         end do
      end do
   end do lp

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)

   ! Add trace of the A matrix
   do iat = 1, mol%nat
      dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
   end do

   if (any(abs(dadr(:, :, :) - numgrad(:, :, :)) > thr_)) then
      call test_failed(error, "Derivative of the A matrix does not match")
      print"(a)", "dadr:"
      print"(3es21.12)", dadr
      print"(a)", "numgrad:"
      print"(3es21.12)", numgrad
      print"(a)", "diff:"
      print"(3es21.12)", dadr - numgrad
   end if

end subroutine test_dadr

subroutine test_dadL(error, mol, model, thr_in)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: ic, jc, iat
   real(wp) :: thr_
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 5.0e-6_wp, unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: xyz(:, :)
   real(wp), allocatable :: qvec(:), numsigma(:, :, :), amatr(:, :), amatl(:, :)
   real(wp) :: eps(3, 3)
   type(cache_container), allocatable :: cache
   allocate (cache)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & amatr(mol%nat + 1, mol%nat + 1), amatl(mol%nat + 1, mol%nat + 1), &
      & dadr(3, mol%nat, mol%nat + 1), dadL(3, 3, mol%nat + 1), atrace(3, mol%nat), &
      & numsigma(3, 3, mol%nat + 1), qvec(mol%nat), xyz(3, mol%nat))

   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%solve(mol, error, cn, qloc, qvec=qvec)
   if (allocated(error)) return

   numsigma = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lp: do ic = 1, 3
      do jc = 1, 3
         amatr(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatr)
         if (allocated(error)) exit lp

         amatl(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         do iat = 1, mol%nat
            ! Numerical sigma of the a matrix
            numsigma(jc, ic, :) = numsigma(jc, ic, :) + &
               & 0.5_wp*qvec(iat)*(amatr(iat, :) - amatl(iat, :))/step
         end do
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)

   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)
   if (allocated(error)) return

   if (any(abs(dadL(:, :, :) - numsigma(:, :, :)) > thr_)) then
      call test_failed(error, "Derivative of the A matrix does not match")
      print"(a)", "dadL:"
      print"(3es21.12)", dadL
      print"(a)", "numsigma:"
      print"(3es21.12)", numsigma
      print"(a)", "diff:"
      print"(3es21.12)", dadL - numsigma
   end if

end subroutine test_dadL

subroutine test_dbdr(error, mol, model, thr_in)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   real(wp) :: thr_
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 5.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dbdr(:, :, :), dbdL(:, :, :)
   real(wp), allocatable :: numgrad(:, :, :), xvecr(:), xvecl(:)
   type(cache_container), allocatable :: cache
   allocate (cache)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & xvecr(mol%nat + 1), xvecl(mol%nat + 1), numgrad(3, mol%nat, mol%nat + 1), &
      & dbdr(3, mol%nat, mol%nat + 1), dbdL(3, 3, mol%nat + 1))

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         ! Right-hand side
         xvecr(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat, :) = 0.5_wp*(xvecr(:) - xvecl(:))/step
      end do
   end do lp

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   call model%get_xvec(mol, cache, xvecl) ! need to call this for xtmp in cache (eeqbc)
   call model%get_xvec_derivs(mol, cache, dbdr, dbdL)

   if (any(abs(dbdr(:, :, :) - numgrad(:, :, :)) > thr_)) then
      call test_failed(error, "Derivative of the b vector does not match")
      print"(a)", "dbdr:"
      print"(3es21.14)", dbdr
      print"(a)", "numgrad:"
      print"(3es21.14)", numgrad
      print"(a)", "diff:"
      print"(3es21.14)", dbdr - numgrad
   end if

end subroutine test_dbdr

subroutine test_dbdL(error, mol, model, thr_in)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic, jc
   real(wp) :: thr_
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 5.0e-6_wp, unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dbdr(:, :, :), dbdL(:, :, :)
   real(wp), allocatable :: numsigma(:, :, :), xvecr(:), xvecl(:)
   real(wp), allocatable :: xyz(:, :)
   real(wp) :: eps(3, 3)
   type(cache_container), allocatable :: cache
   allocate (cache)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & xvecr(mol%nat + 1), xvecl(mol%nat + 1), numsigma(3, 3, mol%nat + 1), &
      & dbdr(3, mol%nat, mol%nat + 1), dbdL(3, 3, mol%nat + 1), xyz(3, mol%nat))

   numsigma = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lp: do ic = 1, 3
      do jc = 1, 3
         ! Right-hand side
         xvecr(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecl)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         do iat = 1, mol%nat
            numsigma(jc, ic, iat) = 0.5_wp*(xvecr(iat) - xvecl(iat))/step
         end do
      end do
   end do lp

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   call model%get_xvec(mol, cache, xvecl) ! need to call this for xtmp in cache (eeqbc)
   call model%get_xvec_derivs(mol, cache, dbdr, dbdL)

   if (any(abs(dbdL(:, :, :) - numsigma(:, :, :)) > thr_)) then
      call test_failed(error, "Derivative of the b vector does not match")
      print"(a)", "dbdL:"
      print"(3es21.14)", dbdL
      print"(a)", "numsigma:"
      print"(3es21.14)", numsigma
      print"(a)", "diff:"
      print"(3es21.14)", dbdL - numsigma
   end if

end subroutine test_dbdL

subroutine gen_test(error, mol, model, qref, eref, efield, thr_in)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Reference charges
   real(wp), intent(in), optional :: qref(:)

   !> Reference energies
   real(wp), intent(in), optional :: eref(:)

   !> Optional external electric field
   real(wp), intent(in), contiguous, optional :: efield(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   real(wp) :: thr_
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: cn(:), qloc(:)
   real(wp), allocatable :: energy(:)
   real(wp), allocatable :: qvec(:)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate (cn(mol%nat), qloc(mol%nat))

   call model%ncoord%get_coordination_number(mol, trans, cn)
   if (allocated(model%ncoord_en)) then
      call model%local_charge(mol, trans, qloc)
   end if

   if (present(eref)) then
      allocate (energy(mol%nat))
      energy(:) = 0.0_wp
   end if
   if (present(qref)) then
      allocate (qvec(mol%nat))
   end if

   call model%solve(mol, error, cn, qloc, energy=energy, qvec=qvec, efield=efield)
   if (allocated(error)) return

   if (present(qref)) then
      if (any(abs(qvec - qref) > thr_)) then
         call test_failed(error, "Partial charges do not match")
         print"(a)", "Charges:"
         print"(3es21.14)", qvec
         print"(a)", "diff:"
         print"(3es21.14)", qvec - qref
      end if
   end if
   if (allocated(error)) return

   if (present(eref)) then
      if (any(abs(energy - eref) > thr_)) then
         call test_failed(error, "Energies do not match")
         print"(a)", "Energy:"
         print"(3es21.14)", energy
         print"(a)", "diff:"
         print"(3es21.14)", energy - eref
      end if
   end if

end subroutine gen_test

subroutine test_numgrad(error, mol, model, thr_in)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   real(wp) :: thr_
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 5.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: numgrad(:, :)
   real(wp) :: er, el

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), sigma(3, 3), numgrad(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         energy(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, gradient=gradient, sigma=sigma)
   if (allocated(error)) return

   if (any(abs(gradient(:, :) - numgrad(:, :)) > thr_)) then
      call test_failed(error, "Derivative of energy does not match")
      print"(a)", "Energy gradient:"
      print"(3es21.14)", gradient
      print"(a)", "numgrad:"
      print"(3es21.14)", numgrad
      print"(a)", "diff:"
      print"(3es21.14)", gradient - numgrad
   end if

end subroutine test_numgrad

subroutine test_numsigma(error, mol, model, thr_in)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: ic, jc
   real(wp) :: thr_
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 5.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :)
   real(wp), allocatable :: xyz(:, :)
   real(wp) :: er, el, eps(3, 3), numsigma(3, 3), sigma(3, 3)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), xyz(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lp: do ic = 1, 3
      do jc = 1, 3
         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         numsigma(jc, ic) = 0.5_wp*(er - el)/step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   energy(:) = 0.0_wp
   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, energy, gradient, sigma)
   if (allocated(error)) return

   if (any(abs(sigma(:, :) - numsigma(:, :)) > thr_)) then
      call test_failed(error, "Derivative of energy does not match")
      print"(a)", "Energy strain:"
      print"(3es21.14)", sigma
      print"(a)", "numsigma:"
      print"(3es21.14)", numsigma
      print"(a)", "diff:"
      print"(3es21.14)", sigma - numsigma
   end if

end subroutine test_numsigma

subroutine test_numdqdr(error, mol, model, thr_in)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   real(wp) :: thr_
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 5.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: ql(:), qr(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & ql(mol%nat), qr(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=qr)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=ql)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp*(qr - ql)/step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdr(:, :, :) - numdr(:, :, :)) > thr_)) then
      call test_failed(error, "Derivative of charges does not match")
      print"(a)", "Charge gradient:"
      print"(3es21.14)", dqdr
      print"(a)", "numgrad:"
      print"(3es21.14)", numdr
      print"(a)", "diff:"
      print"(3es21.14)", dqdr - numdr
   end if

end subroutine test_numdqdr

subroutine test_numdqdL(error, mol, model, thr_in)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: ic, jc
   real(wp) :: thr_
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 5.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: qr(:), ql(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :), numdL(:, :, :)
   real(wp) :: eps(3, 3)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & qr(mol%nat), ql(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & xyz(3, mol%nat), numdL(3, 3, mol%nat))

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=qr)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=ql)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp*(qr - ql)/step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdL(:, :, :) - numdL(:, :, :)) > thr_)) then
      call test_failed(error, "Derivative of charges does not match")
      print"(a)", "Charge gradient:"
      print"(3es21.14)", dqdL
      print"(a)", "numgrad:"
      print"(3es21.14)", numdL
      print"(a)", "diff:"
      print"(3es21.14)", dqdL - numdL
   end if

end subroutine test_numdqdL

subroutine test_eeq_dadr_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model, thr_in=thr1)

end subroutine test_eeq_dadr_mb01

subroutine test_eeq_dadL_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model, thr_in=thr1)

end subroutine test_eeq_dadL_mb01

subroutine test_eeq_dbdr_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model, thr_in=thr1*10)

end subroutine test_eeq_dbdr_mb01

subroutine test_eeq_dbdL_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdL(error, mol, model, thr_in=thr1*10)

end subroutine test_eeq_dbdL_mb01

subroutine test_eeq_q_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & 7.73347759615437E-1_wp, 1.07626897257271E-1_wp, -3.66999554268267E-1_wp, &
      & 4.92833775451616E-2_wp, -1.83332153007808E-1_wp, 2.33302084420314E-1_wp, &
      & 6.61837602813735E-2_wp, -5.43944147972069E-1_wp, -2.70264297953247E-1_wp, &
      & 2.66618970100409E-1_wp, 2.62725030332215E-1_wp, -7.15315061940473E-2_wp, &
      &-3.73300836681036E-1_wp, 3.84585142200261E-2_wp, -5.05851076468890E-1_wp, &
      & 5.17677178773158E-1_wp]

   real(wp), allocatable :: qvec(:)

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)
   if (allocated(error)) return

   ! Check wrapper functions
   allocate (qvec(mol%nat), source=0.0_wp)
   call get_charges(model, mol, error, qvec)
   if (allocated(error)) return

   if (any(abs(qvec - ref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print"(a)", "Charges:"
      print"(3es21.14)", qvec
      print"(a)", "diff:"
      print"(3es21.14)", qvec - ref
   end if
   if (allocated(error)) return

   qvec = 0.0_wp
   call get_eeq_charges(mol, error, qvec)
   if (allocated(error)) return

   if (any(abs(qvec - ref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print"(a)", "Charges:"
      print"(3es21.14)", qvec
      print"(a)", "diff:"
      print"(3es21.14)", qvec - ref
   end if
   if (allocated(error)) return

end subroutine test_eeq_q_mb01

subroutine test_eeq_q_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & 7.38394752482521E-2_wp, -1.68354859084778E-1_wp, -3.47642846358022E-1_wp, &
      &-7.05489251302223E-1_wp, 7.73548241620680E-1_wp, 2.30207580650128E-1_wp, &
      & 1.02748505731185E-1_wp, 9.47818154871089E-2_wp, 2.44259536057649E-2_wp, &
      & 2.34984928231320E-1_wp, -3.17839956573785E-1_wp, 6.67112952465234E-1_wp, &
      &-4.78119957747208E-1_wp, 6.57536208287042E-2_wp, 1.08259091466373E-1_wp, &
      &-3.58215294268738E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeq_q_mb02

subroutine test_eeq_q_efield_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & -4.64683071797936E-2_wp, -5.88362906740707E-1_wp, -2.43193732528000E-1_wp, &
      &  3.51507795107742E-1_wp,  1.07024670739935E+0_wp, -1.33168464005268E+0_wp, &
      & -3.16175726821207E-1_wp,  1.00361461209188E-1_wp, -4.01096542813355E-1_wp, &
      &  1.92924990021589E-1_wp,  5.17680704422172E-1_wp, -7.08449332513941E-1_wp, &
      & -6.26659781146154E-1_wp,  2.23247991084403E-1_wp,  1.68482074733023E+0_wp, &
      &  1.21300573221157E-1_wp]

   ! Electric field
   real(wp), parameter :: efield(3) = [&
      & 0.2000000000000_wp, 0.00000000000000_wp, 0.00000000000000_wp]
      
   call get_structure(mol, "MB16-43", "03")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref, efield=efield)

end subroutine test_eeq_q_efield_mb03

subroutine test_eeq_q_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(17) = [&
      & 1.86904766283711E-02_wp, 2.89972818160259E-01_wp, 3.59298070941105E-02_wp, &
      &-4.61256458126589E-02_wp, -7.02605348653647E-02_wp, -7.42052215689073E-02_wp, &
      &-8.21938718945845E-02_wp, 1.64953118841151E-01_wp, 2.10381640633390E-01_wp, &
      &-6.65485355096282E-02_wp, -2.59873890255450E-01_wp, 1.33839147940414E-01_wp, &
      & 7.20768968601809E-02_wp, -3.36652347675997E-03_wp, -1.14546280789657E-01_wp, &
      &-8.55922398441004E-02_wp, -1.23131162140762E-01_wp]

   call get_structure(mol, "f-block", "Fr_to_Lr")

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeq_q_actinides

subroutine test_eeq_e_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      &-2.18698345033562E-1_wp, -1.04793885931268E+0_wp, 4.78963353574572E-2_wp, &
      & 5.76566377591676E-1_wp, 7.37187470977927E-1_wp, 8.06020047053305E-2_wp, &
      &-4.19837955782898E-1_wp, 5.49627510550566E-2_wp, 8.01486728591565E-2_wp, &
      & 1.00618944521776E-1_wp, -6.61715169034150E-1_wp, -3.60531647289563E-1_wp, &
      &-4.87729666337974E-1_wp, 2.48257554279938E-1_wp, 6.96027176590956E-1_wp, &
      & 4.31679925875087E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeq_e_mb03

subroutine test_eeq_e_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & 1.13974214746111E-1_wp, -4.41735365367827E-1_wp, 8.99036489938394E-2_wp, &
      &-2.97539904703271E-1_wp, 8.05174117097006E-3_wp, 1.31105783760276E-1_wp, &
      & 1.54594451996644E-1_wp, 1.19929653841255E-1_wp, 1.26056586309101E-1_wp, &
      & 1.78439005754586E-1_wp, -1.98703462666082E-1_wp, 4.19630120027785E-1_wp, &
      & 7.05569220334930E-2_wp, -4.50925107441869E-1_wp, 1.39289602382354E-1_wp, &
      &-2.67853086061429E-1_wp]

   call get_structure(mol, "MB16-43", "04")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeq_e_mb04

subroutine test_eeq_g_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model, thr_in=thr1*10)

end subroutine test_eeq_g_mb05

subroutine test_eeq_g_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "06")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model, thr_in=thr1*10)

end subroutine test_eeq_g_mb06

subroutine test_eeq_s_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "07")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model, thr_in=thr1*10)

end subroutine test_eeq_s_mb07

subroutine test_eeq_s_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "08")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model, thr_in=thr1*10)

end subroutine test_eeq_s_mb08

subroutine test_eeq_dqdr_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "09")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model, thr_in=thr1*100)

end subroutine test_eeq_dqdr_mb09

subroutine test_eeq_dqdr_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "10")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model, thr_in=thr1*10)

end subroutine test_eeq_dqdr_mb10

subroutine test_eeq_dqdL_mb11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "11")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model, thr_in=thr1*10)

end subroutine test_eeq_dqdL_mb11

subroutine test_eeq_dqdL_mb12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "12")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model, thr_in=thr1*10)

end subroutine test_eeq_dqdL_mb12

subroutine test_g_h2plus(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   integer, parameter :: nat = 2
   real(wp), parameter :: charge = 1.0_wp
   integer, parameter :: num(nat) = [1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & +0.00000000000000_wp, +0.00000000000000_wp, +0.00000000000000_wp, &
      & +1.00000000000000_wp, +0.00000000000000_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model, thr_in=thr1*100)

end subroutine test_g_h2plus

subroutine test_eeq_dadr_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
      & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
      & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
      & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model, thr_in=thr1)

end subroutine test_eeq_dadr_znooh

subroutine test_eeq_dbdr_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
      & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
      & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
      & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model, thr_in=thr1)

end subroutine test_eeq_dbdr_znooh

subroutine test_g_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
      & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
      & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
      & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model, thr_in=thr1*10)

end subroutine test_g_znooh

subroutine test_dqdr_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
      & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
      & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
      & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model, thr_in=thr1*10)

end subroutine test_dqdr_znooh

subroutine test_eeqbc_dadr_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_dadr_mb01

subroutine test_eeqbc_dadL_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_dadL_mb01

subroutine test_eeqbc_dbdr_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_dbdr_mb01

subroutine test_eeqbc_dbdL_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdL(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_dbdL_mb01

subroutine test_eeqbc_dadr_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_dadr_mb05

subroutine test_eeqbc_dadL_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_dadL_mb05

subroutine test_eeqbc_dbdr_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_dbdr_mb05

subroutine test_eeqbc_q_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      &  6.32665177635486E-1_wp, -7.79000021394748E-3_wp, -7.28560267968198E-1_wp, &
      & -4.10042373033769E-2_wp, -4.60026363641469E-1_wp,  2.20004260276381E-1_wp, &
      & -4.51954034871547E-2_wp, -7.57357565689225E-1_wp, -5.27004217539908E-1_wp, &
      &  3.09566865491112E-1_wp,  2.27494681048596E-1_wp, -3.45550645258975E-1_wp, &
      & -3.99729823174914E-2_wp,  8.48742713551116E-1_wp, -4.49898375805073E-1_wp, &
      &  1.16388636122213E+0_wp]

   real(wp), allocatable :: qvec(:)

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

   ! Check wrapper functions
   allocate (qvec(mol%nat), source=0.0_wp)
   call get_charges(model, mol, error, qvec)
   if (allocated(error)) return

   if (any(abs(qvec - ref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print"(a)", "Charges:"
      print"(3es21.14)", qvec
      print"(a)", "diff:"
      print"(3es21.14)", qvec - ref
   end if
   if (allocated(error)) return

   qvec = 0.0_wp
   call get_eeqbc_charges(mol, error, qvec)
   if (allocated(error)) return

   if (any(abs(qvec - ref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print"(a)", "Charges:"
      print"(3es21.14)", qvec
      print"(a)", "diff:"
      print"(3es21.14)", qvec - ref
   end if
   if (allocated(error)) return

end subroutine test_eeqbc_q_mb01

subroutine test_eeqbc_q_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & -1.68706023409942E-2_wp, -4.31585754567721E-1_wp, -5.63457264808067E-3_wp, &
      & -6.88598679593678E-1_wp,  8.95507266279121E-1_wp,  2.70560711967719E-1_wp, &
      &  1.06704514101413E-2_wp, -1.87359068857452E-3_wp,  2.75916162155502E-1_wp, &
      &  2.68351655726932E-1_wp,  4.83228968763888E-3_wp,  5.84779537866739E-1_wp, &
      & -6.15231036823021E-1_wp,  1.47042131652275E-2_wp, -1.44439879995398E-2_wp, &
      & -5.51084063597411E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeqbc_q_mb02

subroutine test_eeqbc_q_efield_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & -9.77656344417366E-2_wp, -7.23208027030783E-1_wp, -2.16042671980982E-1_wp, &
      &  2.15695910621611E-1_wp,  1.33432687215416E+0_wp, -7.54000025942417E-1_wp, &
      & -4.31828632235353E-1_wp,  1.06805985793153E-1_wp, -3.31034861875571E-1_wp, &
      &  2.11773287076685E-1_wp,  4.58249586177574E-1_wp, -7.03325216892509E-1_wp, &
      & -1.02212584040853E+0_wp,  2.67199973674787E-1_wp,  1.57262943024298E+0_wp, &
      &  1.12649865066930E-1_wp]

   ! Electric field
   real(wp), parameter :: efield(3) = [&
      & 0.2000000000000_wp, 0.00000000000000_wp, 0.00000000000000_wp]
      
   call get_structure(mol, "MB16-43", "03")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref, efield=efield)

end subroutine test_eeqbc_q_efield_mb03

subroutine test_eeqbc_q_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(17) = [&
      &  4.80669255571553E-2_wp, -2.23681443476404E-1_wp,  2.95879870118081E-1_wp, &
      &  1.40528551895345E-1_wp,  1.54496956558730E-1_wp, -1.31765078691778E-1_wp, &
      &  6.43024168002473E-2_wp, -3.08519563494370E-1_wp, -2.76459245598841E-1_wp, &
      &  1.78293689441428E-1_wp, -2.11018657500951E-1_wp, -1.03628773361279E-1_wp, &
      & -1.71248308078648E-1_wp,  2.54400229067594E-1_wp, -5.83023049918706E-2_wp, &
      &  2.01328580342047E-1_wp,  1.47326155413513E-1_wp]

   call get_structure(mol, "f-block", "Fr_to_Lr")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeqbc_q_actinides

subroutine test_eeqbc_e_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      &-6.85383909094970E-02_wp,-2.19074341545678E+00_wp,-4.74601764501974E-03_wp, &
      &-7.87937923890688E-01_wp,-3.03527925157768E+00_wp,-2.82220139737265E-01_wp, &
      &-5.40497002004504E-01_wp,-5.41240617940733E-03_wp,-8.05028961511096E-03_wp, &
      &-1.75953554576722E-02_wp,-1.37304680376243E+00_wp,-3.53988490287759E-01_wp, &
      &-9.94227733791934E-01_wp,-1.58035207823740E-01_wp,-1.51308961030806E+00_wp, &
      &-4.01533041771599E-03_wp]

   call get_structure(mol, "MB16-43", "03")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeqbc_e_mb03

subroutine test_eeqbc_e_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      &-3.11006035589390E-03_wp,-5.28599433581519E-02_wp,-9.35709515842972E-05_wp, &
      &-6.16264514697450E-01_wp,-1.87462481346083E+00_wp,-2.21662455395278E-03_wp, &
      &-2.51039203171614E-02_wp,-8.93995363870797E-03_wp,-1.62270805124076E-04_wp, &
      &-9.27331161591396E-04_wp,-4.55713551016470E-01_wp,-1.49805355046533E+00_wp, &
      &-1.06116580759546E-02_wp,-1.45309833473282E+00_wp,-2.30657216723778E-02_wp, &
      & 1.53766210354731E-06_wp]

   call get_structure(mol, "MB16-43", "04")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeqbc_e_mb04

subroutine test_eeqbc_g_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_g_mb05

subroutine test_eeqbc_g_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "06")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_g_mb06

subroutine test_eeqbc_s_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "07")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_s_mb07

subroutine test_eeqbc_s_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "08")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_s_mb08

subroutine test_eeqbc_dqdr_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "09")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model, thr_in=thr1*10)

end subroutine test_eeqbc_dqdr_mb09

subroutine test_eeqbc_dqdr_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "10")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model, thr_in=thr1*100)

end subroutine test_eeqbc_dqdr_mb10

subroutine test_eeqbc_dqdL_mb11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "11")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model, thr_in=thr1*10)

end subroutine test_eeqbc_dqdL_mb11

subroutine test_eeqbc_dqdL_mb12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "12")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model, thr_in=thr1*10)

end subroutine test_eeqbc_dqdL_mb12

end module test_model
