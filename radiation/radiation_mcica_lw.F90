! radiation_mcica_lw.F90 - Monte-Carlo Independent Column Approximation longtwave solver
!
! (C) Copyright 2015- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Modifications
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-07-12  R. Hogan  Call fast adding method if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_mcica_lw

    use omptimer, only           : omptimer_mark
 use parkind1,              only : jprb
 use radiation_pdf_sampler, only : pdf_sampler_type
 use random_numbers_mix,    only : randomnumberstream, &
      initialize_random_numbers, uniform_distribution


  public

contains

  subroutine generate_column_exp_ran_lr(ng, nlev, ig, random_stream, pdf_sampler, &
    &  frac, pair_cloud_cover, &
    &  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
    &  itrigger, iend, od_scaling)

 use parkind1,              only : jprb
 use radiation_pdf_sampler, only : pdf_sampler_type
 use random_numbers_mix,    only : randomnumberstream, &
      initialize_random_numbers, uniform_distribution


 implicit none

 ! Number of g points / columns, and number of current column
 integer, intent(in) :: ng, ig

 ! Number of levels
 integer, intent(in) :: nlev

 ! Stream for producing random numbers
 type(randomnumberstream), intent(inout) :: random_stream

 ! Object for sampling from a lognormal or gamma distribution
 type(pdf_sampler_type), intent(in) :: pdf_sampler

 ! Cloud fraction, cumulative cloud cover and fractional standard
 ! deviation in each layer
 real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std

 ! Cloud cover of a pair of layers, and amount by which cloud at
 ! next level increases total cloud cover as seen from above
 real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

 ! Overlap parameter of inhomogeneities
 real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom

 ! Top of highest cloudy layer (in this subcolumn) and base of
 ! lowest
 integer, intent(in) :: itrigger, iend

 ! Optical depth scaling to output
 real(jprb), intent(inout), dimension(nlev,ng) :: od_scaling

 ! Height indices
 integer :: jlev, jcloud

 ! Number of contiguous cloudy layers for which to compute optical
 ! depth scaling
 integer :: n_layers_to_scale

 integer :: iy

 ! Is it time to fill the od_scaling variable?
 logical :: do_fill_od_scaling

 real(jprb) :: rand_cloud(nlev)
 real(jprb) :: rand_inhom1(nlev), rand_inhom2(nlev)

 real(jprb) :: omphook_handle

 call omptimer_mark('generate_column_exp_ran_lr',0,omphook_handle)


 ! So far our vertically contiguous cloud contains only one layer
 n_layers_to_scale = 1
 iy = 0

 ! Locate the clouds below this layer: first generate some more
 ! random numbers
 call uniform_distribution(rand_cloud(1:(iend+1-itrigger)),random_stream)

 ! Loop from the layer below the local cloud top down to the
 ! bottom-most cloudy layer
 do jlev = itrigger+1,iend+1
   do_fill_od_scaling = .false.
   if (jlev <= iend) then
     iy = iy+1
     if (n_layers_to_scale > 0) then
       ! There is a cloud above, in which case the probability
       ! of cloud in the layer below is as follows
       if (rand_cloud(iy)*frac(jlev-1) &
            &  < frac(jlev) + frac(jlev-1) - pair_cloud_cover(jlev-1)) then
         ! Add another cloudy layer
         n_layers_to_scale = n_layers_to_scale + 1
       else 
         ! Reached the end of a contiguous set of cloudy layers and
         ! will compute the optical depth scaling immediately.
         do_fill_od_scaling = .true.
       end if
     else
       ! There is clear-sky above, in which case the
       ! probability of cloud in the layer below is as follows
       if (rand_cloud(iy)*(cum_cloud_cover(jlev-1) - frac(jlev-1)) &
            &  < pair_cloud_cover(jlev-1) - overhang(jlev-1) - frac(jlev-1)) then
         ! A new cloud top
         n_layers_to_scale = 1
       end if
     end if
   else
     ! We are at the bottom of the cloudy layers in the model,
     ! so in a moment need to populate the od_scaling array
     do_fill_od_scaling = .true.
   end if

   if (do_fill_od_scaling) then
     ! We have a contiguous range of layers for which we
     ! compute the od_scaling elements using some random
     ! numbers
     call uniform_distribution(rand_inhom1(1:n_layers_to_scale),random_stream)
     call uniform_distribution(rand_inhom2(1:n_layers_to_scale),random_stream)

     ! Loop through the sequence of cloudy layers
     do jcloud = 2,n_layers_to_scale
       ! Use second random number, and inhomogeneity overlap
       ! parameter, to decide whether the first random number
       ! should be repeated (corresponding to maximum overlap)
       ! or not (corresponding to random overlap)
       if (rand_inhom2(jcloud) &
            &  < overlap_param_inhom(jlev-n_layers_to_scale+jcloud-2)) then
         rand_inhom1(jcloud) = rand_inhom1(jcloud-1)
       end if
     end do
     
     ! Sample from a lognormal or gamma distribution to obtain
     ! the optical depth scalings
     call pdf_sampler%sample(fractional_std(jlev-n_layers_to_scale:jlev-1), &
          & rand_inhom1(1:n_layers_to_scale), od_scaling(jlev-n_layers_to_scale:jlev-1,ig))

     n_layers_to_scale = 0
   end if
       
 end do

 call omptimer_mark('generate_column_exp_ran_lr',1,omphook_handle)

end subroutine generate_column_exp_ran_lr



  subroutine generate_column_exp_exp_lr(ng, nlev, ig, random_stream, pdf_sampler, &
    &  frac, pair_cloud_cover, &
    &  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
    &  itrigger, iend, od_scaling)

 use parkind1,              only : jprb
 use radiation_pdf_sampler, only : pdf_sampler_type
 use random_numbers_mix,    only : randomnumberstream, &
      initialize_random_numbers, uniform_distribution

 implicit none

 ! Number of g points / columns, and number of current column
 integer, intent(in) :: ng, ig

 ! Number of levels
 integer, intent(in) :: nlev

 ! Stream for producing random numbers
 type(randomnumberstream), intent(inout) :: random_stream

 ! Object for sampling from a lognormal or gamma distribution
 type(pdf_sampler_type), intent(in) :: pdf_sampler

 ! Cloud fraction, cumulative cloud cover and fractional standard
 ! deviation in each layer
 real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std

 ! Cloud cover of a pair of layers, and amount by which cloud at
 ! next level increases total cloud cover as seen from above
 real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

 ! Overlap parameter of inhomogeneities
 real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom

 ! Top of highest cloudy layer (in this subcolumn) and base of
 ! lowest
 integer, intent(in) :: itrigger, iend

 ! Optical depth scaling to output
 real(jprb), intent(inout), dimension(nlev,ng) :: od_scaling

 ! Height indices
 integer :: jlev, jcloud

 integer :: iy

 real(jprb) :: rand_cloud(nlev)
 real(jprb) :: rand_inhom1(nlev), rand_inhom2(nlev)

 ! For each column analysed, this vector locates the clouds. It is
 ! only actually used for Exp-Exp overlap
 logical :: is_cloudy(nlev)

 ! Number of contiguous cloudy layers for which to compute optical
 ! depth scaling
 integer :: n_layers_to_scale

 real(jprb) :: omphook_handle

 call omptimer_mark('generate_column_exp_exp_lr',0,omphook_handle)

 iy = 0

 is_cloudy = .false.
 is_cloudy(itrigger) = .true.

 ! Locate the clouds below this layer: first generate some more
 ! random numbers
 call uniform_distribution(rand_cloud(1:(iend+1-itrigger)),random_stream)

 ! Loop from the layer below the local cloud top down to the
 ! bottom-most cloudy layer
 do jlev = itrigger+1,iend
   iy = iy+1
   if (is_cloudy(jlev-1)) then
     ! There is a cloud above, in which case the probability
     ! of cloud in the layer below is as follows
     if (rand_cloud(iy)*frac(jlev-1) &
          &  < frac(jlev) + frac(jlev-1) - pair_cloud_cover(jlev-1)) then
       ! Add another cloudy layer
       is_cloudy(jlev) = .true.
     end if
   else
     ! There is clear-sky above, in which case the
     ! probability of cloud in the layer below is as follows
     if (rand_cloud(iy)*(cum_cloud_cover(jlev-1) - frac(jlev-1)) &
          &  < pair_cloud_cover(jlev-1) - overhang(jlev-1) - frac(jlev-1)) then
         ! A new cloud top
       is_cloudy(jlev) = .true.
     end if
   end if
 end do

 ! We have a contiguous range of layers for which we compute the
 ! od_scaling elements using some random numbers

 ! In the Exp-Exp overlap scheme we do all layers at once
 n_layers_to_scale = iend+1 - itrigger
     
 call uniform_distribution(rand_inhom1(1:n_layers_to_scale),random_stream)
 call uniform_distribution(rand_inhom2(1:n_layers_to_scale),random_stream)
     
 ! Loop through the sequence of cloudy layers
 do jcloud = 2,n_layers_to_scale
   ! Use second random number, and inhomogeneity overlap
   ! parameter, to decide whether the first random number
   ! should be repeated (corresponding to maximum overlap)
   ! or not (corresponding to random overlap)
   if (rand_inhom2(jcloud) &
        &  < overlap_param_inhom(iend-n_layers_to_scale+jcloud-1)) then
     rand_inhom1(jcloud) = rand_inhom1(jcloud-1)
   end if
 end do
     
 ! Sample from a lognormal or gamma distribution to obtain the
 ! optical depth scalings, calling the faster masked version and
 ! assuming values outside the range itrigger:iend are already zero
 call pdf_sampler%masked_sample(n_layers_to_scale, &
      &  fractional_std(itrigger:iend), &
      &  rand_inhom1(1:n_layers_to_scale), od_scaling(itrigger:iend,ig), &
      &  is_cloudy(itrigger:iend))
     
 call omptimer_mark('generate_column_exp_exp_lr',1,omphook_handle)

end subroutine generate_column_exp_exp_lr



  !---------------------------------------------------------------------
  ! Longwave Monte Carlo Independent Column Approximation
  ! (McICA). This implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. This method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. The cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
subroutine solver_mcica_lw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od_in, ssa_in, g_in, od_cloud_in, ssa_cloud_in, g_cloud_in, planck_hl_in, &
       &  emission_in, albedo_in, &
       &  flux)

  use parkind1, only           : jprb,jprd
  use yomhook,  only           : lhook, dr_hook
  use radiation_io,   only           : nulerr, radiation_abort
  use radiation_config, only         : config_type
  use radiation_single_level, only   : single_level_type
  use radiation_cloud, only          : cloud_type
  use radiation_cloud_cover, only : IOverlapExponential
  use radiation_flux, only           : flux_type
  use radiation_two_stream, only     : calc_two_stream_gammas_lw, &
      &                               calc_two_stream_gammas_lw_lr, &
      &                               calc_two_stream_gammas_lw_cond_lr, &
      &                               calc_reflectance_transmittance_lw, &
      &                               calc_reflectance_transmittance_lw_lr, &
      &                               calc_reflectance_transmittance_lw_cond_lr, &
      &                               calc_no_scattering_transmittance_lw, &
      &                               calc_no_scattering_transmittance_lw_lr, &
      &                               calc_no_scattering_transmittance_lw_cond_lr, &
#ifdef FAST_EXPONENTIAL
         &                               exp_fast, &
#endif
         &                               LwDiffusivity

#ifndef FAST_EXPONENTIAL
#define exp_fast exp
#endif

  use radiation_adding_ica_lw, only  : adding_ica_lw_lr, adding_ica_lw_cond_lr, fast_adding_ica_lw, &
      &                               fast_adding_ica_lw_lr, calc_fluxes_no_scattering_lw, &
      &                               calc_fluxes_no_scattering_lw_lr, calc_fluxes_no_scattering_lw_cond_lr
  use radiation_lw_derivatives, only : calc_lw_derivatives_ica, calc_lw_derivatives_ica_lr, &
  &                                    modify_lw_derivatives_ica, modify_lw_derivatives_ica_lr2
  use radiation_cloud_generator, only: cloud_generator_lr !, generate_column_exp_ran_lr, generate_column_exp_exp_lr

  use random_numbers_mix, only : randomnumberstream

  implicit none

  ! Inputs
  integer, intent(in) :: nlev               ! number of model levels
  integer, intent(in) :: istartcol, iendcol ! range of columns to process
  type(config_type),        intent(in) :: config
  type(single_level_type),  intent(in) :: single_level
  type(cloud_type),         intent(in) :: cloud

  ! Gas and aerosol optical depth, single-scattering albedo and
  ! asymmetry factor at each longwave g-point
  real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
      &  od_in
  real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
      &  ssa_in, g_in

  ! Cloud and precipitation optical depth, single-scattering albedo and
  ! asymmetry factor in each longwave band
  real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
      &  od_cloud_in
  real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
      &  nlev,istartcol:iendcol) :: ssa_cloud_in, g_cloud_in

  ! Planck function at each half-level and the surface
  real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
      &  planck_hl_in

  ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
  ! surface at each longwave g-point
  real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission_in, albedo_in

  ! Output
  type(flux_type), intent(inout):: flux

  ! Local variables

  ! Diffuse reflectance and transmittance for each layer in clear
  ! and all skies
  real(jprb), dimension(istartcol:iendcol, nlev) :: ref_clear, reflectance
  ! cos: ng can not be demoted because of reductions
  real(jprb), dimension(istartcol:iendcol,nlev) :: transmittance

  real(jprb), dimension(istartcol:iendcol,nlev) :: trans_clear

  ! Emission by a layer into the upwelling or downwelling diffuse
  ! streams, in clear and all skies
  ! cos: ng can not be demoted because of reductions
  real(jprb), dimension(istartcol:iendcol,nlev) :: source_up_clear, source_up, &
                                              source_dn_clear, source_dn

  ! Fluxes per g point
  ! cos: ng can not be demoted because of reductions
  real(jprb), dimension(istartcol:iendcol,nlev+1) :: flux_dn, flux_dn_clear, flux_up, flux_up_clear, &
  &                   flux_up_clear_sum, flux_up_sum, flux_dn_sum, flux_dn_clear_sum
  real(jprb), dimension(istartcol:iendcol,nlev) :: flux_up_mul_trans_clear_sum, flux_up_mul_trans_sum
  real(jprb) :: flux_up_mul_trans_clear_prod, flux_up_mul_trans_prod

  ! Combined gas+aerosol+cloud optical depth, single scattering
  ! albedo and asymmetry factor
  real(jprb), dimension(istartcol:iendcol) :: ssa_total, od_total, g_total

  ! Combined scattering optical depth
  real(jprb) :: scat_od, scat_od_total

  ! Two-stream coefficients
  ! cos: could be demoted to scalar if calc_two_stream_gammas_lw_lr and 
  ! calc_reflectance_transmittance_lw_lr are fused in the same jloop 
  real(jprb), dimension(istartcol:iendcol) :: gamma1, gamma2 

  ! Optical depth scaling from the cloud generator, zero indicating
  ! clear skies
  ! cos; (ng) can be demoted, but required moving the last ng loop of cloud_generator out and fuse it
  real(jprb), dimension(istartcol:iendcol,nlev,config%n_g_lw) :: od_scaling

  ! Modified optical depth after McICA scaling to represent cloud
  ! inhomogeneity
  ! cos: temporarily added jcol & nlev for loop splitting
  real(jprb), dimension(istartcol:iendcol) :: od_cloud_new

  ! Total cloud cover output from the cloud generator
  real(jprb), dimension(istartcol:iendcol) :: total_cloud_cover

  ! Identify clear-sky layers
  logical :: is_clear_sky_layer(istartcol:iendcol,nlev)

  ! Index of the highest cloudy layer
  ! cos : temporarily added jcol
  integer :: i_cloud_top(istartcol:iendcol)

  real(jprb) :: factor  
  
  ! Number of g points
  integer :: ng

  ! Loop indices for level, column and g point
  integer :: jlev, jcol, jg

  real(jprb) :: hook_handle, omphook_solver_mcica_lw, omphook_calc_two_stream_gammas_lw, &
  &  omphook_adding_ica_lw, omphook_calc_no_scattering_transmittance_lw, &
  & omphook_calc_fluxes_no_scattering_lw, omphook_cloud_generator, &
  & omphook_set_scat_od, omphook_calc_two_stream_gammas_lw_b, &
  & omphook_calc_reflectance_transmittance_lw, &
  & omphook_calc_no_scattering_transmittance_lw_b, &
  & omphook_adding_ica_lw_b, omphook_fast_adding_ica_lw, &
  & omphook_calc_fluxes_no_scattering_lw_b, omphook_calc_lw_derivatives_ica, &
  & omphook_modify_lw_derivatives_ica, omphook_generate_column_exp_exp, &
  & omphook_calc_lw_derivatives_ica_lr, omphook_marker1, omphook_marker2, omphook_marker3           

  real(jprb), dimension(istartcol:iendcol, nlev, config%n_g_lw) :: &
  &  od

  real(jprb), dimension(istartcol:iendcol, nlev, config%n_g_lw_if_scattering) :: &
  &  ssa, g

  ! Cloud and precipitation optical depth, single-scattering albedo and
  ! asymmetry factor in each longwave band
  real(jprb), dimension(istartcol:iendcol,nlev,config%n_bands_lw)   :: &
      &  od_cloud
  real(jprb), dimension(istartcol:iendcol,nlev,config%n_bands_lw_if_scattering) :: ssa_cloud, g_cloud

  ! Planck function at each half-level and the surface
  real(jprb), dimension(istartcol:iendcol,nlev+1,config%n_g_lw) :: &
      &  planck_hl

  ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
  ! surface at each longwave g-point
  real(jprb), dimension(istartcol:iendcol, config%n_g_lw) :: emission, albedo

  ! Scaled random number for finding cloud
  real(jprb) :: trigger
  integer :: itrigger

  ! Uniform deviates between 0 and 1
  ! cos: original (ng). Future demote to jcol only
  real(jprb) :: rand_top(istartcol:iendcol,config%n_g_lw)

  ! First and last cloudy layers
  ! cos: origina (scalar). Need to remain like that
  integer :: ibegin(istartcol:iendcol), iend(istartcol:iendcol)

  ! Cloud cover of a pair of layers, and amount by which cloud at
  ! next level increases total cloud cover as seen from above
  ! cos: original (nlev+1). Future can not be demoted since we move jcol innermost
  ! we need to retain the independent jcol calculations
  real(jprb), dimension(istartcol:iendcol, nlev-1) :: pair_cloud_cover, overhang

  ! Cumulative cloud cover from TOA to the base of each layer
  ! cos: original (lev). Future can not be demoted since we move jcol innermost
  ! we need to retain the independent jcol calculations
  real(jprb) :: cum_cloud_cover(istartcol:iendcol,nlev)

      ! Overlap parameter of inhomogeneities
  !cos: original (nlev). Future can not be demoted
  real(jprb) :: overlap_param_inhom(istartcol:iendcol,nlev-1)

  ! Seed for random number generator and stream for producing random
  ! numbers
  !cos: oritinal (scalar)
  type(randomnumberstream) :: random_stream(istartcol:iendcol)


  real(jprd) :: coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot                          
  real(jprd) :: k_exponent, reftrans_factor
  real(jprd) :: exponential  ! = exp(-k_exponent*od)
  real(jprd) :: exponential2 ! = exp(-2*k_exponent*od)

  ! cos: temporaries for fast_adding_ica_lw_lr
  ! ---------------------------------------------
  ! Albedo of the entire earth/atmosphere system below each half
  ! level
  real(jprb), dimension(istartcol:iendcol,nlev+1) :: albedo_tmp

  ! Upwelling radiation at each half-level due to emission below
  ! that half-level (W m-2)
  real(jprb), dimension(istartcol:iendcol,nlev+1) :: source

  ! Equal to 1/(1-albedo*reflectance)
  real(jprb), dimension(istartcol:iendcol,nlev)   :: inv_denominator
  ! ---------------------------------------------

  integer, dimension(iendcol-istartcol+1) :: cloud_cover_idx
  integer :: cloud_cover_idx_length, idx

  integer, dimension(iendcol-istartcol+1,nlev) :: cloud_cover_fraction
  integer, dimension(nlev) :: cloud_cover_fraction_length

  ! Height indices
  integer :: jcloud

  integer :: iy

  real(jprb) :: rand_cloud(nlev)
  real(jprb) :: rand_inhom1(nlev), rand_inhom2(nlev)

  ! For each column analysed, this vector locates the clouds. It is
  ! only actually used for Exp-Exp overlap
  logical :: is_cloudy(nlev)

  ! Number of contiguous cloudy layers for which to compute optical
  ! depth scaling
  integer :: n_layers_to_scale

  ! Is it time to fill the od_scaling variable?
  logical :: do_fill_od_scaling

  ng = config%n_g_lw

  write(*,*) "DOMAIN SIZES :", ng, nlev, istartcol,iendcol

  do jcol = istartcol,iendcol
    do jlev=1,nlev
      do jg = 1,config%n_g_lw
            od(jcol,jlev,jg) = od_in(jg,jlev,jcol)
      enddo
    enddo
  enddo

  do jcol = istartcol,iendcol
    do jlev=1,nlev
      do jg = 1,config%n_g_lw_if_scattering
            ssa(jcol,jlev,jg) = ssa_in(jg,jlev,jcol)
            g(jcol,jlev,jg) = g_in(jg,jlev,jcol)
      enddo
    enddo
  enddo

  do jcol = istartcol,iendcol
    do jlev=1,nlev
      do jg = 1,config%n_bands_lw
            od_cloud(jcol,jlev,jg) = od_cloud_in(jg,jlev,jcol)
      enddo
    enddo
  enddo

  do jcol = istartcol,iendcol
    do jlev=1,nlev
      do jg = 1,config%n_bands_lw_if_scattering
            ssa_cloud(jcol,jlev,jg) = ssa_cloud_in(jg,jlev,jcol)
            g_cloud(jcol,jlev,jg) = g_cloud_in(jg,jlev,jcol)
      enddo
    enddo
  enddo

  do jcol = istartcol,iendcol
    do jlev=1,nlev+1
      do jg = 1,config%n_g_lw
            planck_hl(jcol,jlev,jg) = planck_hl_in(jg,jlev,jcol)
      enddo
    enddo
  enddo

  do jcol = istartcol,iendcol
    do jg = 1,config%n_g_lw
      emission(jcol,jg) = emission_in(jg,jcol)
      albedo(jcol,jg) = albedo_in(jg,jcol)
    enddo
  enddo

  if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',0,hook_handle)
call omptimer_mark('radiation_mcica_lw:solver_mcica_lw',0, &
     &  omphook_solver_mcica_lw)


  if (.not. config%do_clear) then
    write(nulerr,'(a)') '*** Error: longwave McICA requires clear-sky calculation to be performed'
    call radiation_abort()      
  end if

!call omptimer_mark('cloud_generator',0, &
!&   omphook_cloud_generator)

  ! Do cloudy-sky calculation; add a prime number to the seed in
  ! the longwave
  call cloud_generator_lr(ng, istartcol, iendcol, nlev, config%i_overlap_scheme, &
        &  single_level%iseed, &
        &  config%cloud_fraction_threshold, &
        &  cloud%fraction, cloud%overlap_param, &
        &  config%cloud_inhom_decorr_scaling, &
        &  total_cloud_cover, rand_top, &
        &  pair_cloud_cover, overhang, cum_cloud_cover, overlap_param_inhom, &
        &  random_stream, ibegin, iend, &
        &  is_beta_overlap=config%use_beta_overlap)
      

!call omptimer_mark('cloud_generator',1, &
!&   omphook_cloud_generator)

  
  do jcol = istartcol, iendcol
      ! Store total cloud cover
      flux%cloud_cover_lw(jcol) = total_cloud_cover(jcol)      
  enddo

  cloud_cover_idx = 0

  cloud_cover_idx_length = 0
  do jcol = istartcol,iendcol
      if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
      cloud_cover_idx_length = cloud_cover_idx_length + 1
      cloud_cover_idx(cloud_cover_idx_length) = jcol
      endif
  enddo

  cloud_cover_fraction = 0
  do jlev = 1,nlev
      cloud_cover_fraction_length(jlev) = 0
      do jcol = istartcol,iendcol
        if ( (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) .and. &
  &          (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold)) then
          cloud_cover_fraction_length(jlev) = cloud_cover_fraction_length(jlev)+1
          cloud_cover_fraction(cloud_cover_fraction_length(jlev), jlev) = jcol
        endif
      enddo
  enddo

  do idx = 1, cloud_cover_idx_length
      jcol = cloud_cover_idx(idx)
      ! Reset optical depth scaling to clear skies
      od_scaling(jcol,:,:) = 0.0_jprb

      is_clear_sky_layer(jcol,:) = .true.
      i_cloud_top(jcol) = nlev+1
  enddo

  do jlev = 1,nlev
    do idx = 1, cloud_cover_fraction_length(jlev)
      jcol = cloud_cover_fraction(idx,jlev)
      is_clear_sky_layer(jcol,jlev) = .false.
      ! Get index to the first cloudy layer from the top
      if (i_cloud_top(jcol) > jlev) then
        i_cloud_top(jcol) = jlev
      endif
    enddo
  enddo
      
  flux_up_clear_sum(:,:) = 0.0
  flux_dn_clear_sum(:,:) = 0.0

  flux_up_sum(:,:) = 0.0
  flux_dn_sum(:,:) = 0.0

  flux_up_mul_trans_clear_sum(:,:) = 0.0
  flux_up_mul_trans_sum(:,:) = 0.0

  do jg = 1,ng
!    call omptimer_mark('generate_column_exp_exp',0,omphook_generate_column_exp_exp)
!    call omptimer_mark('cloud_generator',0, &
!          &   omphook_cloud_generator)

    do idx = 1, cloud_cover_idx_length
      jcol = cloud_cover_idx(idx)
      ! Loop over ng columns
      ! cos: the random num generation was before out of the innermost loop. 
      ! With the reordering had to be brought inside. We would need to refactor 
      ! the random number generation subroutines to recover performance

      ! Find the cloud top height corresponding to the current
      ! random number, and store in itrigger
      trigger = rand_top(jcol,jg) * total_cloud_cover(jcol)
      jlev = ibegin(jcol)
      do while (trigger > cum_cloud_cover(jcol,jlev) .and. jlev < iend(jcol))
        jlev = jlev + 1
      end do
      itrigger = jlev

      if (config%i_overlap_scheme /= IOverlapExponential) then
     ! cos: inline of generate_column_exp_ran_lr
     !            call generate_column_exp_ran_lr(ng, nlev, jg, random_stream(jcol), config%pdf_sampler, &
     !                &  cloud%fraction(jcol,:), pair_cloud_cover(jcol,:), &
     !                &  cum_cloud_cover(jcol,:), overhang(jcol,:), cloud%fractional_std(jcol,:), overlap_param_inhom(jcol,:), &
     !                &  itrigger, iend(jcol), od_scaling(jcol,:,:))

        ! So far our vertically contiguous cloud contains only one layer
        n_layers_to_scale = 1
        iy = 0

        ! Locate the clouds below this layer: first generate some more
        ! random numbers
        call uniform_distribution(rand_cloud(1:(iend(jcol)+1-itrigger)),random_stream(jcol))

        ! Loop from the layer below the local cloud top down to the
        ! bottom-most cloudy layer
        do jlev = itrigger+1,iend(jcol)+1
          do_fill_od_scaling = .false.
          if (jlev <= iend(jcol)) then
            iy = iy+1
            if (n_layers_to_scale > 0) then
              ! There is a cloud above, in which case the probability
              ! of cloud in the layer below is as follows
              if (rand_cloud(iy)*cloud%fraction(jcol,jlev-1) &
                  &  < cloud%fraction(jcol,jlev) + cloud%fraction(jcol,jlev-1) - pair_cloud_cover(jcol,jlev-1)) then
                  ! Add another cloudy layer
                n_layers_to_scale = n_layers_to_scale + 1
              else
                ! Reached the end of a contiguous set of cloudy layers and
                ! will compute the optical depth scaling immediately.
                do_fill_od_scaling = .true.
              end if
            else
              ! There is clear-sky above, in which case the
              ! probability of cloud in the layer below is as follows
              if (rand_cloud(iy)*(cum_cloud_cover(jcol,jlev-1) - cloud%fraction(jcol,jlev-1)) &
                  &  < pair_cloud_cover(jcol,jlev-1) - overhang(jcol,jlev-1) - cloud%fraction(jcol,jlev-1)) then
                  ! A new cloud top
                n_layers_to_scale = 1
              end if
            end if
          else
            ! We are at the bottom of the cloudy layers in the model,
            ! so in a moment need to populate the od_scaling array
            do_fill_od_scaling = .true.
          end if

          if (do_fill_od_scaling) then
            ! We have a contiguous range of layers for which we
            ! compute the od_scaling elements using some random
            ! numbers
            call uniform_distribution(rand_inhom1(1:n_layers_to_scale),random_stream(jcol))
            call uniform_distribution(rand_inhom2(1:n_layers_to_scale),random_stream(jcol))

            ! Loop through the sequence of cloudy layers
            do jcloud = 2,n_layers_to_scale
              ! Use second random number, and inhomogeneity overlap
              ! parameter, to decide whether the first random number
              ! should be repeated (corresponding to maximum overlap)
              ! or not (corresponding to random overlap)
              if (rand_inhom2(jcloud) &
                  &  < overlap_param_inhom(jcol,jlev-n_layers_to_scale+jcloud-2)) then
                rand_inhom1(jcloud) = rand_inhom1(jcloud-1)
              end if
            end do
            
            ! Sample from a lognormal or gamma distribution to obtain
            ! the optical depth scalings
            call config%pdf_sampler%sample(cloud%fractional_std(jcol,jlev-n_layers_to_scale:jlev-1), &
                  & rand_inhom1(1:n_layers_to_scale), od_scaling(jcol,jlev-n_layers_to_scale:jlev-1,jg))

            n_layers_to_scale = 0
          end if
        end do

        ! cos: END inline of generate_column_exp_ran_lr

      else
        call generate_column_exp_exp_lr(ng, nlev, jg, random_stream(jcol), config%pdf_sampler, &
                         &  cloud%fraction(jcol,:), pair_cloud_cover(jcol,:), &
                         &  cum_cloud_cover(jcol,:), overhang(jcol,:), cloud%fractional_std(jcol,:), overlap_param_inhom(jcol,:), &
                         &  itrigger, iend(jcol), od_scaling(jcol,:,:))
      end if      
    end do
!call omptimer_mark('generate_column_exp_exp',1,omphook_generate_column_exp_exp)
!call omptimer_mark('cloud_generator',1, &
!               &   omphook_cloud_generator)


    ! Clear-sky calculation
    if (config%do_lw_aerosol_scattering) then
      ! Scattering case: first compute clear-sky reflectance,
      ! transmittance etc at each model level
!call omptimer_mark('calc_two_stream_gammas_lw',0, &
!            &  omphook_calc_two_stream_gammas_lw) 
      do jlev = 1,nlev
        ssa_total = ssa(:,jlev,jg)
        g_total   = g(:,jlev,jg)
! cos: inline of calc_two_stream_gammas_lw_lr
!          call calc_two_stream_gammas_lw_lr(istartcol, iendcol, ssa_total(:), g_total, &
!               &  gamma1, gamma2)

        do jcol = istartcol,iendcol
          ! Fu et al. (1997), Eq 2.9 and 2.10:
          !      gamma1(jg) = LwDiffusivity * (1.0_jprb - 0.5_jprb*ssa(jg) &
          !           &                    * (1.0_jprb + g(jg)))
          !      gamma2(jg) = LwDiffusivity * 0.5_jprb * ssa(jg) &
          !           &                    * (1.0_jprb - g(jg))
          ! Reduce number of multiplications
          factor = (LwDiffusivity * 0.5_jprb) * ssa_total(jcol)
          gamma1(jcol) = LwDiffusivity - factor*(1.0_jprb + g_total(jcol))
          gamma2(jcol) = factor * (1.0_jprb - g_total(jcol))
        end do
! cos: END inline of calc_two_stream_gammas_lw_lr

! cos: inline of calc_reflectance_transmittance_lw_lr
!          call calc_reflectance_transmittance_lw_lr(istartcol, iendcol, &
!               &  od(:,jlev,jg), gamma1, gamma2, &
!               &  planck_hl(:,jlev,jg), planck_hl(:,jlev+1,jg), &
!               &  ref_clear(:,jlev), trans_clear(:,jlev), &
!               &  source_up_clear(:,jlev), source_dn_clear(:,jlev))

        do jcol = istartcol,iendcol
          if (od(jcol,jlev,jg) > 1.0e-3_jprd) then
            k_exponent = sqrt(max((gamma1(jcol) - gamma2(jcol)) * (gamma1(jcol) + gamma2(jcol)), &
                  1.E-12_jprd)) ! Eq 18 of Meador & Weaver (1980)
            exponential = exp_fast(-k_exponent*od(jcol,jlev,jg))
            exponential2 = exponential*exponential
            reftrans_factor = 1.0 / (k_exponent + gamma1(jcol) + (k_exponent - gamma1(jcol))*exponential2)
            ! Meador & Weaver (1980) Eq. 25
            ref_clear(jcol,jlev) = gamma2(jcol) * (1.0_jprd - exponential2) * reftrans_factor
            ! Meador & Weaver (1980) Eq. 26
            trans_clear(jcol,jlev) = 2.0_jprd * k_exponent * exponential * reftrans_factor
        
            ! Compute upward and downward emission assuming the Planck
            ! function to vary linearly with optical depth within the layer
            ! (e.g. Wiscombe , JQSRT 1976).

            ! Stackhouse and Stephens (JAS 1991) Eqs 5 & 12
            coeff = (planck_hl(jcol,jlev+1,jg)-planck_hl(jcol,jlev,jg)) / (od(jcol,jlev,jg)*(gamma1(jcol)+gamma2(jcol)))
            coeff_up_top  =  coeff + planck_hl(jcol,jlev,jg)
            coeff_up_bot  =  coeff + planck_hl(jcol,jlev+1,jg)
            coeff_dn_top  = -coeff + planck_hl(jcol,jlev,jg)
            coeff_dn_bot  = -coeff + planck_hl(jcol,jlev+1,jg)
            source_up_clear(jcol,jlev) =  coeff_up_top - ref_clear(jcol,jlev) * coeff_dn_top - trans_clear(jcol,jlev) * coeff_up_bot
            source_dn_clear(jcol,jlev) =  coeff_dn_bot - ref_clear(jcol,jlev) * coeff_up_bot - trans_clear(jcol,jlev) * coeff_dn_top
          else
            k_exponent = sqrt(max((gamma1(jcol) - gamma2(jcol)) * (gamma1(jcol) + gamma2(jcol)), &
                  1.E-12_jprd)) ! Eq 18 of Meador & Weaver (1980)
            ref_clear(jcol,jlev) = gamma2(jcol) * od(jcol,jlev,jg)
            trans_clear(jcol,jlev) = (1.0_jprb - k_exponent*od(jcol,jlev,jg)) / &
            &     (1.0_jprb + od(jcol,jlev,jg)*(gamma1(jcol)-k_exponent))
            source_up_clear(jcol,jlev) = (1.0_jprb - ref_clear(jcol,jlev) - trans_clear(jcol,jlev)) &
                  &       * 0.5 * (planck_hl(jcol,jlev,jg) + planck_hl(jcol,jlev+1,jg))
            source_dn_clear(jcol,jlev) = source_up_clear(jcol,jlev)
          end if
        end do
 ! cos: END inline of calc_reflectance_transmittance_lw_lr
      end do
!call omptimer_mark('calc_two_stream_gammas_lw',1, &
!      &  omphook_calc_two_stream_gammas_lw) 

!call omptimer_mark('adding_ica_lw',0, &
!      &   omphook_adding_ica_lw) 

      ! Then use adding method to compute fluxes
      call adding_ica_lw_lr(istartcol, iendcol, nlev, &
          &  ref_clear, trans_clear(:,:), source_up_clear, source_dn_clear, &
          &  emission(:,jg), albedo(:,jg), &
          &  flux_up_clear(:,:), flux_dn_clear(:,:))

!call omptimer_mark('adding_ica_lw',1, &
!      &   omphook_adding_ica_lw)         

    else

!call omptimer_mark('calc_no_scattering_transmittance_lw',0, &
!          &   omphook_calc_no_scattering_transmittance_lw)

          ! ! Non-scattering case: use simpler functions for
          ! ! transmission and emission
      do jlev = 1,nlev
! cos: inline of calc_no_scattering_transmittance_lw_lr
!          call calc_no_scattering_transmittance_lw_lr(istartcol, iendcol, od(:,jlev,jg), &
!               &  planck_hl(:,jlev,jg), planck_hl(:,jlev+1,jg), &
!               &  trans_clear(:,jlev), source_up_clear(:,jlev), source_dn_clear(:,jlev))
        do jcol = istartcol,iendcol  
          ! Compute upward and downward emission assuming the Planck
          ! function to vary linearly with optical depth within the layer
          ! (e.g. Wiscombe , JQSRT 1976).
          if (od(jcol,jlev,jg) > 1.0e-3) then
            ! Simplified from calc_reflectance_transmittance_lw above
            coeff = LwDiffusivity*od(jcol,jlev,jg)
            trans_clear(jcol,jlev) = exp_fast(-coeff)
            coeff = (planck_hl(jcol,jlev+1,jg)-planck_hl(jcol,jlev,jg)) / coeff
            coeff_up_top  =  coeff + planck_hl(jcol,jlev,jg)
            coeff_up_bot  =  coeff + planck_hl(jcol,jlev+1,jg)
            coeff_dn_top  = -coeff + planck_hl(jcol,jlev,jg)
            coeff_dn_bot  = -coeff + planck_hl(jcol,jlev+1,jg)
            source_up_clear(jcol,jlev) =  coeff_up_top - trans_clear(jcol,jlev) * coeff_up_bot
            source_dn_clear(jcol,jlev) =  coeff_dn_bot - trans_clear(jcol,jlev) * coeff_dn_top
          else
            ! Linear limit at low optical depth
            coeff = LwDiffusivity*od(jcol,jlev,jg)
            trans_clear(jcol,jlev) = 1.0_jprb - coeff
            source_up_clear(jcol,jlev) = coeff * 0.5_jprb * (planck_hl(jcol,jlev,jg)+planck_hl(jcol,jlev+1,jg))
            source_dn_clear(jcol,jlev) = source_up_clear(jcol,jlev)
          end if
        end do
! cos: END inline of calc_no_scattering_transmittance_lw_lr
      end do

!call omptimer_mark('calc_no_scattering_transmittance_lw',1, &
!          &   omphook_calc_no_scattering_transmittance_lw)

!call omptimer_mark('calc_fluxes_no_scattering_lw',0, &
!          &   omphook_calc_fluxes_no_scattering_lw)


! cos: inline of calc_fluxes_no_scattering_lw_lr
!        ! ! ! Simpler down-then-up method to compute fluxes
!        call calc_fluxes_no_scattering_lw_lr(istartcol, iendcol, nlev, &
!             &  trans_clear(:,:), source_up_clear, source_dn_clear, &
!             &  emission(:,jg), albedo(:,jg), &
!             &  flux_up_clear(:,:), flux_dn_clear(:,:))
          ! At top-of-atmosphere there is no diffuse downwelling radiation
      flux_dn_clear(:,1) = 0.0_jprb

      ! Work down through the atmosphere computing the downward fluxes
      ! at each half-level
      do jlev = 1,nlev
        flux_dn_clear(:,jlev+1) = trans_clear(:,jlev)*flux_dn_clear(:,jlev) + source_dn_clear(:,jlev)
      end do

      ! Surface reflection and emission
      flux_up_clear(:,nlev+1) = emission(:,jg) + albedo(:,jg) * flux_dn_clear(:,nlev+1)

      ! Work back up through the atmosphere computing the upward fluxes
      ! at each half-level
      do jlev = nlev,1,-1
        flux_up_clear(:,jlev) = trans_clear(:,jlev)*flux_up_clear(:,jlev+1) + source_up_clear(:,jlev)
      end do

            
!call omptimer_mark('calc_fluxes_no_scattering_lw',1, &
!      &   omphook_calc_fluxes_no_scattering_lw)



      ! Ensure that clear-sky reflectance is zero since it may be
      ! used in cloudy-sky case
      ref_clear = 0.0_jprb
      end if

      do jlev = 1,nlev
        do idx = 1,cloud_cover_fraction_length(jlev)
          jcol = cloud_cover_fraction(idx, jlev)
          od_cloud_new(jcol) = od_scaling(jcol,jlev,jg) &
                &  * od_cloud(jcol,jlev,config%i_band_from_reordered_g_lw(jg))
          od_total(jcol) = od(jcol,jlev,jg) + od_cloud_new(jcol)
          ssa_total(jcol) = 0.0_jprb
          g_total(jcol)   = 0.0_jprb
        enddo

        if (config%do_lw_cloud_scattering) then
          ! Scattering case: calculate reflectance and
          ! transmittance at each model level

!call omptimer_mark('set_scat_od',0, &
!               &   omphook_set_scat_od)

          do idx = 1,cloud_cover_fraction_length(jlev)
            jcol = cloud_cover_fraction(idx, jlev)

            if (config%do_lw_aerosol_scattering) then
              ! In single precision we need to protect against the
              ! case that od_total > 0.0 and ssa_total > 0.0 but
              ! od_total*ssa_total == 0 due to underflow
              scat_od_total = ssa(jcol,jlev,jg)*od(jcol,jlev,jg) &
                  &     + ssa_cloud(jcol,jlev,config%i_band_from_reordered_g_lw(jg)) &
                  &     *  od_cloud_new(jcol)
              if (scat_od_total > 0.0_jprb) then
                g_total(jcol) = (g(jg,jlev,jcol)*ssa(jcol,jlev,jg)*od(jcol,jlev,jg) &
                    &     +   g_cloud(jcol,jlev,config%i_band_from_reordered_g_lw(jg)) &
                    &     * ssa_cloud(jcol,jlev,config%i_band_from_reordered_g_lw(jg)) &
                    &     *  od_cloud_new(jcol)) &
                    &     / scat_od_total
              endif                
              if (od_total(jcol) > 0.0_jprb) then
               ssa_total(jcol) = scat_od_total / od_total(jcol)
              endif
            else
              if (od_total(jcol) > 0.0_jprb) then
                scat_od = ssa_cloud(jcol,jlev,config%i_band_from_reordered_g_lw(jg)) &
                &     * od_cloud_new(jcol)
                ssa_total(jcol) = scat_od / od_total(jcol)
                if (scat_od > 0.0_jprb) then
                  g_total(jcol) = g_cloud(jcol,jlev,config%i_band_from_reordered_g_lw(jg)) &
                      &     * ssa_cloud(jcol,jlev,config%i_band_from_reordered_g_lw(jg)) &
                      &     *  od_cloud_new(jcol) / scat_od
                end if
              end if
                    !end do
            end if
          enddo
      
!call omptimer_mark('set_scat_od',1, &
!               &   omphook_set_scat_od)


!call omptimer_mark('calc_two_stream_gammas_lw_b',0, &
!               &   omphook_calc_two_stream_gammas_lw_b)

               ! Compute cloudy-sky reflectance, transmittance etc at
               ! each model level
               ! cos: inline of calc_two_stream_gammas_lw_cond_lr
               ! call calc_two_stream_gammas_lw_cond_lr(istartcol, iendcol, total_cloud_cover, cloud%fraction(:,jlev), &
               !                  & config%cloud_fraction_threshold, ssa_total, g_total, gamma1, gamma2)

          do idx = 1,cloud_cover_fraction_length(jlev)
            jcol = cloud_cover_fraction(idx, jlev)
            
            ! Fu et al. (1997), Eq 2.9 and 2.10:
            !      gamma1(jg) = LwDiffusivity * (1.0_jprb - 0.5_jprb*ssa(jg) &
            !           &                    * (1.0_jprb + g(jg)))
            !      gamma2(jg) = LwDiffusivity * 0.5_jprb * ssa(jg) &
            !           &                    * (1.0_jprb - g(jg))
            ! Reduce number of multiplications
            factor = (LwDiffusivity * 0.5_jprb) * ssa_total(jcol)
            gamma1(jcol) = LwDiffusivity - factor*(1.0_jprb + g_total(jcol))
            gamma2(jcol) = factor * (1.0_jprb - g_total(jcol))
          end do

!call omptimer_mark('calc_two_stream_gammas_lw_b',1, &
!&   omphook_calc_two_stream_gammas_lw_b)

!call omptimer_mark('calc_reflectance_transmittance_lw',0, &
!&   omphook_calc_reflectance_transmittance_lw)


! cos: inline of calc_reflectance_transmittance_lw_cond_lr
!          call calc_reflectance_transmittance_lw_cond_lr(istartcol, iendcol, &
!                  &  total_cloud_cover, cloud%fraction(:,jlev), config%cloud_fraction_threshold, &
!                  & od_total, gamma1, gamma2, &
!                  &  planck_hl(:,jlev,jg), planck_hl(:,jlev+1,jg), &
!                  &  reflectance(:,jlev), transmittance(:,jlev), &
!                  & source_up(:,jlev), source_dn(:,jlev))

          do idx = 1,cloud_cover_fraction_length(jlev)
            jcol = cloud_cover_fraction(idx, jlev)
                
            if (od_total(jcol) > 1.0e-3_jprd) then
              k_exponent = sqrt(max((gamma1(jcol) - gamma2(jcol)) * (gamma1(jcol) + gamma2(jcol)), &
                    1.E-12_jprd)) ! Eq 18 of Meador & Weaver (1980)
              exponential = exp_fast(-k_exponent*od_total(jcol))
              exponential2 = exponential*exponential
              reftrans_factor = 1.0 / (k_exponent + gamma1(jcol) + (k_exponent - gamma1(jcol))*exponential2)
              ! Meador & Weaver (1980) Eq. 25
              reflectance(jcol,jlev) = gamma2(jcol) * (1.0_jprd - exponential2) * reftrans_factor
              ! Meador & Weaver (1980) Eq. 26
              transmittance(jcol,jlev) = 2.0_jprd * k_exponent * exponential * reftrans_factor
              
              ! Compute upward and downward emission assuming the Planck
              ! function to vary linearly with optical depth within the layer
              ! (e.g. Wiscombe , JQSRT 1976).

              ! Stackhouse and Stephens (JAS 1991) Eqs 5 & 12
              coeff = (planck_hl(jcol,jlev+1,jg)-planck_hl(jcol,jlev,jg)) / (od_total(jcol)*(gamma1(jcol)+gamma2(jcol)))
              coeff_up_top  =  coeff + planck_hl(jcol,jlev,jg)
              coeff_up_bot  =  coeff + planck_hl(jcol,jlev+1,jg)
              coeff_dn_top  = -coeff + planck_hl(jcol,jlev,jg)
              coeff_dn_bot  = -coeff + planck_hl(jcol,jlev+1,jg)
              source_up(jcol,jlev) =  coeff_up_top - reflectance(jcol,jlev) * coeff_dn_top - transmittance(jcol,jlev) * coeff_up_bot
              source_dn(jcol,jlev) =  coeff_dn_bot - reflectance(jcol,jlev) * coeff_up_bot - transmittance(jcol,jlev) * coeff_dn_top
            else
              k_exponent = sqrt(max((gamma1(jcol) - gamma2(jcol)) * (gamma1(jcol) + gamma2(jcol)), &
                    1.E-12_jprd)) ! Eq 18 of Meador & Weaver (1980)
              reflectance(jcol,jlev) = gamma2(jcol) * od_total(jcol)
              transmittance(jcol,jlev) = (1.0_jprb - k_exponent*od_total(jcol)) / & 
              &   (1.0_jprb + od_total(jcol)*(gamma1(jcol)-k_exponent))
              source_up(jcol,jlev) = (1.0_jprb - reflectance(jcol,jlev) - transmittance(jcol,jlev)) &
                    &       * 0.5 * (planck_hl(jcol,jlev,jg) + planck_hl(jcol,jlev+1,jg))
              source_dn(jcol,jlev) = source_up(jcol,jlev)
            end if
          end do
 ! cos: END inline of calc_reflectance_transmittance_lw_cond_lr

!call omptimer_mark('calc_reflectance_transmittance_lw',1, &
!&   omphook_calc_reflectance_transmittance_lw)


        else
!call omptimer_mark('calc_no_scattering_transmittance_lw_b',0, &
!&   omphook_calc_no_scattering_transmittance_lw_b)


          ! No-scattering case: use simpler functions for
          ! transmission and emission
          call calc_no_scattering_transmittance_lw_cond_lr(istartcol, iendcol, &
                  & cloud_cover_fraction(:,jlev), cloud_cover_fraction_length(jlev), &
                  & od_total, planck_hl(:,jlev,jg), planck_hl(:,jlev+1,jg), &
                  &  transmittance(:,jlev), source_up(:,jlev), source_dn(:,jlev))

!call omptimer_mark('calc_no_scattering_transmittance_lw_b',1, &
!&   omphook_calc_no_scattering_transmittance_lw_b)


        end if

        do jcol = istartcol,iendcol

          if ((total_cloud_cover(jcol) >= config%cloud_fraction_threshold) .and. &
&           (cloud%fraction(jcol,jlev) < config%cloud_fraction_threshold)) then

            ! Clear-sky layer: copy over clear-sky values
            reflectance(jcol,jlev) = ref_clear(jcol,jlev)
            transmittance(jcol,jlev) = trans_clear(jcol,jlev)
            source_up(jcol,jlev) = source_up_clear(jcol,jlev)
            source_dn(jcol,jlev) = source_dn_clear(jcol,jlev)
          endif
        enddo
      end do
      if (config%do_lw_aerosol_scattering) then
        ! Use adding method to compute fluxes for an overcast sky,
        ! allowing for scattering in all layers
!call omptimer_mark('adding_ica_lw_b',0, &
!&   omphook_adding_ica_lw_b)

        call adding_ica_lw_cond_lr(istartcol, iendcol, nlev, cloud_cover_idx, cloud_cover_idx_length, &
&          reflectance, transmittance(:,:), source_up, &
&          source_dn, emission(:,jg), albedo(:,jg), flux_up(:,:), flux_dn(:,:))

!call omptimer_mark('adding_ica_lw_b',1, &
!&   omphook_adding_ica_lw_b)

      else if (config%do_lw_cloud_scattering) then
        ! Use adding method to compute fluxes but optimize for the
        ! presence of clear-sky layers
!                      call fast_adding_ica_lw(ng, nlev, &
! &               reflectance(:,:,jcol), transmittance(:,:,jcol), source_up(:,:,jcol), &
!                 & source_dn(:,:,jcol), emission(:,jcol), albedo(:,jcol), is_clear_sky_layer(:,jcol), i_cloud_top(jcol), &
!                 & flux_dn_clear(:,:,jcol), flux_up(:,:,jcol), flux_dn(:,:,jcol))

!call omptimer_mark('fast_adding_ica_lw',0, &
!&   omphook_fast_adding_ica_lw)

! cos: inline of fast_adding_ica_lw_lr
!          call fast_adding_ica_lw_lr(istartcol,iendcol, nlev, total_cloud_cover, config%cloud_fraction_threshold, &
!&               reflectance, transmittance(:,:), source_up, &
!              & source_dn, emission(:,jg), albedo(:,jg), is_clear_sky_layer(:,:), i_cloud_top, &
!              & flux_dn_clear(:,:), flux_up(:,:), flux_dn(:,:))

        do idx = 1,cloud_cover_idx_length
          jcol = cloud_cover_idx(idx)
        ! Copy over downwelling fluxes above cloud from clear sky
          flux_dn(jcol,1:i_cloud_top(jcol)) = flux_dn_clear(jcol,1:i_cloud_top(jcol))

          albedo_tmp(jcol,nlev+1) = albedo(jcol,jg)

          ! At the surface, the source is thermal emission
          source(jcol,nlev+1) = emission(jcol,jg)    
        enddo

      ! Work back up through the atmosphere and compute the albedo of
      ! the entire earth/atmosphere system below that half-level, and
      ! also the "source", which is the upwelling flux due to emission
      ! below that level
      ! do jlev = nlev,i_cloud_top,-1
        do jlev = nlev,1,-1
          do idx = 1,cloud_cover_idx_length
            jcol = cloud_cover_idx(idx)

            if(jlev >= i_cloud_top(jcol)) then
              if (is_clear_sky_layer(jcol,jlev)) then
              ! ! Reflectance of this layer is zero, simplifying the expression

                albedo_tmp(jcol,jlev) = transmittance(jcol,jlev)*transmittance(jcol,jlev)*albedo_tmp(jcol,jlev+1)
                source(jcol,jlev) = source_up(jcol,jlev) &
                      &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
                      &                    + albedo_tmp(jcol,jlev+1)*source_dn(jcol,jlev))
              else
                ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
                inv_denominator(jcol,jlev) = 1.0_jprb &
                      &  / (1.0_jprb-albedo_tmp(jcol,jlev+1)*reflectance(jcol,jlev))
                ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
                albedo_tmp(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev)*transmittance(jcol,jlev) &
                      &  * albedo_tmp(jcol,jlev+1) * inv_denominator(jcol,jlev)
                ! Shonk & Hogan (2008) Eq 11:
                source(jcol,jlev) = source_up(jcol,jlev) &
                      &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
                      &                    + albedo_tmp(jcol,jlev+1)*source_dn(jcol,jlev)) &
                      &                   * inv_denominator(jcol,jlev)
              endif
            endif
        enddo
      enddo
 
      do idx = 1,cloud_cover_idx_length
        jcol = cloud_cover_idx(idx)
        ! Compute the fluxes above the highest cloud
        flux_up(jcol,i_cloud_top(jcol)) = source(jcol,i_cloud_top(jcol)) &
            &                 + albedo_tmp(jcol,i_cloud_top(jcol))*flux_dn(jcol,i_cloud_top(jcol))
      enddo
      !do jlev = i_cloud_top(jcol)-1,1,-1
      ! cos: would it make sense to revert the jcol and jlev loops here
      ! in this pattern to avoid the jcol conditional exit? performance wise need to test?
      do jlev = nlev-1,1,-1
        do idx = 1,cloud_cover_idx_length
          jcol = cloud_cover_idx(idx)
          if(jlev < i_cloud_top(jcol)) then
            flux_up(jcol,jlev) = transmittance(jcol,jlev)*flux_up(jcol,jlev+1) + source_up(jcol,jlev)
          endif
        end do
      enddo

      ! Work back down through the atmosphere from cloud top computing
      ! the fluxes at each half-level

      !do jlev = i_cloud_top,nlev
      do jlev = 1,nlev
        do idx = 1,cloud_cover_idx_length
          jcol = cloud_cover_idx(idx)
          if (jlev >= i_cloud_top(jcol)) then

            if (is_clear_sky_layer(jcol,jlev)) then
              flux_dn(jcol,jlev+1) = transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
                    &               + source_dn(jcol,jlev)
              flux_up(jcol,jlev+1) = albedo_tmp(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
                    &               + source(jcol,jlev+1)
            else
              ! Shonk & Hogan (2008) Eq 14 (after simplification):
              flux_dn(jcol,jlev+1) &
                    &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
                    &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
                    &     + source_dn(jcol,jlev)) * inv_denominator(jcol,jlev)
              ! Shonk & Hogan (2008) Eq 12:
              flux_up(jcol,jlev+1) = albedo_tmp(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
                    &               + source(jcol,jlev+1)
            end if
          endif
        enddo
      enddo

!call omptimer_mark('fast_adding_ica_lw',1, &
!&   omphook_fast_adding_ica_lw)


    else
        ! ! Simpler down-then-up method to compute fluxes
        ! call calc_fluxes_no_scattering_lw(ng, nlev, &
        !      &  transmittance(:,:,jcol), source_up(:,:,jcol), source_dn(:,:,jcol), emission(:,jcol), albedo(:,jcol), &
        !      &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
                        ! Simpler down-then-up method to compute fluxes
!call omptimer_mark('calc_fluxes_no_scattering_lw_b',0, &
!&   omphook_calc_fluxes_no_scattering_lw_b)


      call calc_fluxes_no_scattering_lw_cond_lr(istartcol,iendcol, nlev, cloud_cover_idx, cloud_cover_idx_length, &
        &  transmittance(:,:), source_up, source_dn, emission(:,jg), albedo(:,jg), &
        &  flux_up(:,:), flux_dn(:,:))
!call omptimer_mark('calc_fluxes_no_scattering_lw_b',1, &
!&   omphook_calc_fluxes_no_scattering_lw_b)


    end if

!call omptimer_mark('marker3',0, &
!      &   omphook_marker3)
      
    flux_up_clear_sum(:,:) = flux_up_clear_sum(:,:) + flux_up_clear(:,:)
    flux_dn_clear_sum(:,:) = flux_dn_clear_sum(:,:) + flux_dn_clear(:,:)

    do jcol = istartcol,iendcol
      flux_up_mul_trans_clear_prod  = flux_up_clear(jcol,nlev+1)

      do jlev = nlev,1,-1
        flux_up_mul_trans_clear_prod = flux_up_mul_trans_clear_prod * trans_clear(jcol,jlev)
        flux_up_mul_trans_clear_sum(jcol,jlev) = flux_up_mul_trans_clear_sum(jcol,jlev) + flux_up_mul_trans_clear_prod
      enddo
    enddo

    do idx = 1,cloud_cover_idx_length
      jcol = cloud_cover_idx(idx)
      flux_up_sum(jcol,:) = flux_up_sum(jcol,:) + flux_up(jcol,:)
      flux_dn_sum(jcol,:) = flux_dn_sum(jcol,:) + flux_dn(jcol,:)

      flux_up_mul_trans_prod  = flux_up(jcol,nlev+1)
      do jlev = nlev,1,-1
        flux_up_mul_trans_prod = flux_up_mul_trans_prod * transmittance(jcol,jlev)

        flux_up_mul_trans_sum(jcol,jlev) = flux_up_mul_trans_sum(jcol,jlev) + flux_up_mul_trans_prod
      enddo
    enddo

    do jcol = istartcol,iendcol
      flux%lw_dn_surf_clear_g(jg,jcol) = flux_dn_clear(jcol,nlev+1)
    enddo

    do idx = 1,cloud_cover_idx_length
      jcol = cloud_cover_idx(idx)
            ! Store surface spectral downwelling fluxes
      flux%lw_dn_surf_g(jg,jcol) = total_cloud_cover(jcol)*flux_dn(jcol,nlev+1) &
            &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_dn_surf_clear_g(jg,jcol)
    enddo
!call omptimer_mark('marker3',1, &
!      &   omphook_marker3)
      
  enddo

    ! cos: here there are reductions on ng. Therefore we need to break the ng loop,
    ! and the storages can only be ng independent if the accumulation is performed in previous loops

!call omptimer_mark('marker2',0, &
!    &   omphook_marker2)
    
      ! cos: array syntax once data layouts are compatible
  do jcol = istartcol,iendcol

    ! Sum over g-points to compute broadband fluxes
    flux%lw_up_clear(jcol,:) = flux_up_clear_sum(jcol,:) !
    flux%lw_dn_clear(jcol,:) = flux_dn_clear_sum(jcol,:)
  enddo

!call omptimer_mark('marker2',1, omphook_marker2)
    
  do idx = 1,cloud_cover_idx_length

!call omptimer_mark('marker1',0, omphook_marker1)
    
    jcol = cloud_cover_idx(idx)
      
    ! Store overcast broadband fluxes
    flux%lw_up(jcol,:) = flux_up_sum(jcol,:)
    flux%lw_dn(jcol,:) = flux_dn_sum(jcol,:) 

    ! Cloudy flux profiles currently assume completely overcast
    ! skies; perform weighted average with clear-sky profile
    flux%lw_up(jcol,:) =  total_cloud_cover(jcol) *flux%lw_up(jcol,:) &
          &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_up_clear(jcol,:)
    flux%lw_dn(jcol,:) =  total_cloud_cover(jcol) *flux%lw_dn(jcol,:) &
          &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_dn_clear(jcol,:)

!call omptimer_mark('marker1',1, omphook_marker1)

      ! Compute the longwave derivatives needed by Hogan and Bozzo
      ! (2015) approximate radiation update scheme
    if (config%do_lw_derivatives) then
!call omptimer_mark('calc_lw_derivatives_ica',0, &
!&   omphook_calc_lw_derivatives_ica)


      call calc_lw_derivatives_ica_lr(ng, nlev, jcol, flux_up_sum(jcol,nlev+1), &
        &    flux_up_mul_trans_sum(jcol,:), flux%lw_derivatives)
!call omptimer_mark('calc_lw_derivatives_ica',1, &
!&   omphook_calc_lw_derivatives_ica)



      if (total_cloud_cover(jcol) < 1.0_jprb - config%cloud_fraction_threshold) then
        ! Modify the existing derivative with the contribution from the clear sky
!call omptimer_mark('modify_lw_derivatives_ica',0, &
!&   omphook_modify_lw_derivatives_ica)

        call modify_lw_derivatives_ica_lr2(ng, nlev, jcol, &
&             flux_up_clear_sum(jcol,nlev+1), flux_up_mul_trans_clear_sum(jcol,:), &
&             1.0_jprb-total_cloud_cover(jcol), flux%lw_derivatives)
!call omptimer_mark('modify_lw_derivatives_ica',1, &
!&   omphook_modify_lw_derivatives_ica)


      end if
    end if
  enddo
!call omptimer_mark('calc_lw_derivatives_ica_lr',0, &
!    &  omphook_calc_lw_derivatives_ica_lr) 

  do jcol = istartcol,iendcol
    if (total_cloud_cover(jcol) < config%cloud_fraction_threshold) then
      ! No cloud in profile and clear-sky fluxes already
      ! calculated: copy them over
      flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
      flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
      flux%lw_dn_surf_g(:,jcol) = flux%lw_dn_surf_clear_g(:,jcol)
      if (config%do_lw_derivatives) then
        call calc_lw_derivatives_ica_lr(ng, nlev, jcol, flux_up_clear_sum(jcol,nlev+1), &
        &    flux_up_mul_trans_clear_sum(jcol,:), flux%lw_derivatives)

      end if
    end if ! Cloud is present in profile
  end do

!call omptimer_mark('calc_lw_derivatives_ica_lr',1, &
!&  omphook_calc_lw_derivatives_ica_lr) 

call omptimer_mark('radiation_mcica_lw:solver_mcica_lw',1,omphook_solver_mcica_lw)

if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',1,hook_handle)
    
end subroutine solver_mcica_lw

end module radiation_mcica_lw
