! radiation_mcica_lw_acc.F90 - Monte-Carlo Independent Column Approximation longtwave solver
!
! Copyright (C) 2015-2017 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-07-12  R. Hogan  Call fast adding method if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_mcica_lw_acc

  public

contains

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
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream_acc, only : calc_two_stream_gammas_lw, &
         &                               calc_reflectance_transmittance_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw_acc, only  : adding_ica_lw, fast_adding_ica_lw, &
         &                               calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives_acc, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
    use radiation_cloud_generator, only: cloud_generator

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
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
         &  od_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: ref_clear, trans_clear, reflectance, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: source_up_clear, source_dn_clear, source_up, source_dn

    ! Fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up, flux_dn
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up_clear, flux_dn_clear

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: od_total, ssa_total, g_total

    ! Two-stream coefficients
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: gamma1, gamma2

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: od_scaling

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: od_cloud_new

    ! Total cloud cover output from the cloud generator
    real(jprb), dimension(istartcol:iendcol) :: total_cloud_cover

    ! Identify clear-sky layers
    logical :: is_clear_sky_layer(nlev, istartcol:iendcol)

    ! Index of the highest cloudy layer
    integer :: i_cloud_top(istartcol:iendcol)

    ! logical multi-purpose 2d (nlev, jcol) mask
    logical :: mask_2d(nlev, istartcol:iendcol)
    ! logical multi-purpose 1d (jcol) mask
    logical :: mask_1d(istartcol:iendcol)

    ! Number of g points
    integer :: ng

    ! Loop indices for level and column
    integer :: jlev, jcol, jg

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_lw_acc:solver_mcica_lw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: longwave McICA requires clear-sky calculation to be performed'
      call radiation_abort()      
    end if

    ng = config%n_g_lw

    if (config%do_lw_aerosol_scattering) then
      ssa_total = ssa
      g_total   = g
      mask_2d = .TRUE.
      mask_1d = .TRUE.
      call calc_two_stream_gammas_lw(ng, nlev, istartcol, iendcol, mask_2d, ssa_total, g_total, &
      &  gamma1, gamma2)
      call calc_reflectance_transmittance_lw(ng, nlev, istartcol, iendcol, mask_2d, &
           &  od, gamma1, gamma2, &
           &  planck_hl, &
           &  ref_clear, trans_clear, source_up_clear, source_dn_clear)
      ! Then use adding method to compute fluxes
      call adding_ica_lw(ng, nlev, istartcol, iendcol, mask_1d, &
           &  ref_clear, trans_clear, source_up_clear, source_dn_clear, &
           &  emission, albedo, &
           &  flux_up_clear, flux_dn_clear)

    else
      mask_2d = .TRUE.
      mask_1d = .TRUE.
      call calc_no_scattering_transmittance_lw(ng, nlev, istartcol, iendcol, mask_2d, od, &
      &  planck_hl, &
      &  trans_clear, source_up_clear, source_dn_clear)
      ! Loop through columns
      ! Simpler down-then-up method to compute fluxes
      call calc_fluxes_no_scattering_lw(ng, nlev, istartcol, iendcol, mask_1d, &
      &  trans_clear, source_up_clear, source_dn_clear, &
      &  emission, albedo, &
      &  flux_up_clear, flux_dn_clear)
      ! Ensure that clear-sky reflectance is zero since it may be
      ! used in cloudy-sky case
      ref_clear = 0.0_jprb
    endif

    ! Loop through columns
    do jcol = istartcol,iendcol
      ! Sum over g-points to compute broadband fluxes
      flux%lw_up_clear(jcol,:) = sum(flux_up_clear(:,:,jcol),1)
      flux%lw_dn_clear(jcol,:) = sum(flux_dn_clear(:,:,jcol),1)
      ! Store surface spectral downwelling fluxes
      flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear(:,nlev+1,jcol)
    end do

    ! Loop through columns
    do jcol = istartcol,iendcol
      ! Do cloudy-sky calculation; add a prime number to the seed in
      ! the longwave
      call cloud_generator(ng, nlev, config%i_overlap_scheme, &
           &  single_level%iseed(jcol) + 997, &
           &  config%cloud_fraction_threshold, &
           &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
           &  config%cloud_inhom_decorr_scaling, cloud%fractional_std(jcol,:), &
           &  config%pdf_sampler, od_scaling(:,:,jcol), total_cloud_cover(jcol), &
           &  is_beta_overlap=config%use_beta_overlap)
      
      ! Store total cloud cover
      flux%cloud_cover_lw(jcol) = total_cloud_cover(jcol)
    end do
      
    is_clear_sky_layer(:,:) = .true.
    i_cloud_top(:) = nlev+1
    ! Loop through columns
    do jlev = 1,nlev
      do jcol = istartcol,iendcol
        if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
          ! Compute combined gas+aerosol+cloud optical properties
          if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
            is_clear_sky_layer(jlev,jcol) = .false.
            ! Get index to the first cloudy layer from the top
            if (i_cloud_top(jcol) > jlev) then
              i_cloud_top(jcol) = jlev
            end if
          end if
        end if
      end do
    end do
    ! Loop through columns
    do jlev = 1,nlev
      do jcol = istartcol,iendcol
        do jg = 1, ng
          if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
            ! Compute combined gas+aerosol+cloud optical properties
            if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
              od_cloud_new(jg,jlev,jcol) = od_scaling(jg,jlev,jcol) &
                 &  * od_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol)
              od_total(jg,jlev,jcol) = od(jg,jlev,jcol) + od_cloud_new(jg,jlev,jcol)
            end if
          end if
        end do
      end do
    end do
    ! Loop through columns
    if (config%do_lw_cloud_scattering) then
      ! Scattering case: calculate reflectance and
      ! transmittance at each model level
      if (config%do_lw_aerosol_scattering) then
        do jlev = 1,nlev
          do jcol = istartcol,iendcol
            do jg = 1, ng
              if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
                ! Compute combined gas+aerosol+cloud optical properties
                if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
                  if (od_total(jg,jlev,jcol) > 0.0_jprb) then
                    ssa_total(jg,jlev,jcol) = (ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                       &     + ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                       &     *  od_cloud_new(jg,jlev,jcol)) & 
                       &     / od_total(jg,jlev,jcol)
                  else
                    ssa_total(jg,jlev,jcol) = 0.0_jprb
                  endif
                  if (ssa_total(jg,jlev,jcol) > 0.0_jprb) then
                    g_total(jg,jlev,jcol) = (g(jg,jlev,jcol)*ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                       &     +   g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                       &     *  od_cloud_new(jg,jlev,jcol)) &
                       &     / (ssa_total(jg,jlev,jcol)*od_total(jg,jlev,jcol))
                  else
                    g_total(jg,jlev,jcol) = 0.0_jprb
                  end if
                end if
              end if
            end do
          end do
        end do
      else
        do jlev = 1,nlev
          do jcol = istartcol,iendcol
            do jg = 1, ng
              if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
                ! Compute combined gas+aerosol+cloud optical properties
                if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
                  if (od_total(jg,jlev,jcol) > 0.0_jprb) then
                    ssa_total(jg,jlev,jcol) = ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                       &     * od_cloud_new(jg,jlev,jcol) / od_total(jg,jlev,jcol)
                  else
                    ssa_total(jg,jlev,jcol) = 0.0_jprb
                  endif
                  if (ssa_total(jg,jlev,jcol) > 0.0_jprb) then
                    g_total(jg,jlev,jcol) = g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                       &     *  od_cloud_new(jg,jlev,jcol) / (ssa_total(jg,jlev,jcol)*od_total(jg,jlev,jcol))
                  else
                    g_total(jg,jlev,jcol) = 0.0_jprb
                  endif
                end if
              end if
            end do
          end do
        end do
      end if
      ! Loop through columns
      mask_2d = .FALSE.
      do jlev = 1,nlev
        do jcol = istartcol,iendcol
          if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
            ! Compute combined gas+aerosol+cloud optical properties
            if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
              mask_2d(jlev, jcol) = .TRUE.
            end if
          end if
        end do
      end do
      call calc_two_stream_gammas_lw(ng, nlev, istartcol, iendcol, mask_2d, ssa_total, g_total, &
      &  gamma1, gamma2)
      call calc_reflectance_transmittance_lw(ng, nlev, istartcol, iendcol, mask_2d, &
           &  od_total, gamma1, gamma2, &
           &  planck_hl, &
           &  reflectance, transmittance, source_up, source_dn)
    else
      ! No-scattering case: use simpler functions for
      ! transmission and emission
      call calc_no_scattering_transmittance_lw(ng, nlev, istartcol, iendcol, mask_2d, od_total, &
           &  planck_hl, &
           &  transmittance, source_up, source_dn)
    end if
    ! Loop through columns
    do jlev = 1,nlev
      do jcol = istartcol,iendcol
        if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
          if (cloud%fraction(jcol,jlev) < config%cloud_fraction_threshold) then
            ! Clear-sky layer: copy over clear-sky values
            reflectance(:,jlev,jcol) = ref_clear(:,jlev,jcol)
            transmittance(:,jlev,jcol) = trans_clear(:,jlev,jcol)
            source_up(:,jlev,jcol) = source_up_clear(:,jlev,jcol)
            source_dn(:,jlev,jcol) = source_dn_clear(:,jlev,jcol)
          end if
        end if
      end do
    end do
    mask_1d = .FALSE.
    ! Loop through columns
    do jcol = istartcol,iendcol
      if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
        mask_1d(jcol) = .TRUE.
      end if
    end do
    if (config%do_lw_aerosol_scattering) then
      ! Use adding method to compute fluxes for an overcast sky,
      ! allowing for scattering in all layers
      call adding_ica_lw(ng, nlev, istartcol, iendcol, mask_1d, &
           &  reflectance, transmittance, source_up, source_dn, &
           &  emission, albedo, &
           &  flux_up, flux_dn)
    else if (config%do_lw_cloud_scattering) then
      ! Use adding method to compute fluxes but optimize for the
      ! presence of clear-sky layers
      ! call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
      !      &  emission(:,jcol), albedo(:,jcol), &
      !      &  flux_up, flux_dn)
      call fast_adding_ica_lw(ng, nlev, istartcol, iendcol, mask_1d, &
           &  reflectance, transmittance, source_up, source_dn, &
           &  emission, albedo, &
           &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
           &  flux_up, flux_dn)
    else
      ! Loop through columns
      ! Simpler down-then-up method to compute fluxes
      call calc_fluxes_no_scattering_lw(ng, nlev, istartcol, iendcol, mask_1d, &
           &  transmittance, source_up, source_dn, emission, albedo, &
           &  flux_up, flux_dn)
    end if
    ! Loop through columns
    do jcol = istartcol,iendcol
      if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
        ! Store overcast broadband fluxes
        flux%lw_up(jcol,:) = sum(flux_up(:,:,jcol),1)
        flux%lw_dn(jcol,:) = sum(flux_dn(:,:,jcol),1)

        ! Cloudy flux profiles currently assume completely overcast
        ! skies; perform weighted average with clear-sky profile
        flux%lw_up(jcol,:) =  total_cloud_cover(jcol) *flux%lw_up(jcol,:) &
             &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) =  total_cloud_cover(jcol) *flux%lw_dn(jcol,:) &
             &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_dn_clear(jcol,:)
        ! Store surface spectral downwelling fluxes
        flux%lw_dn_surf_g(:,jcol) = total_cloud_cover(jcol)*flux_dn(:,nlev+1,jcol) &
             &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_dn_surf_clear_g(:,jcol)

        ! Compute the longwave derivatives needed by Hogan and Bozzo
        ! (2015) approximate radiation update scheme
      end if
    end do
    if (config%do_lw_derivatives) then
      call calc_lw_derivatives_ica(ng, nlev, istartcol, iendcol, mask_1d, transmittance(:,:,:), flux_up(:,nlev+1,:), &
           &                       flux%lw_derivatives)
      ! Loop through columns
      mask_1d = .FALSE.
      do jcol = istartcol,iendcol
       if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
         if (total_cloud_cover(jcol) < 1.0_jprb - config%cloud_fraction_threshold) then
           mask_1d(jcol) = .TRUE.
         end if
       end if
      end do
      ! Modify the existing derivative with the contribution from the clear sky
      call modify_lw_derivatives_ica(ng, nlev, istartcol, iendcol, mask_1d, trans_clear(:,:,:), flux_up_clear(:,nlev+1,:), &
           &                         1.0_jprb-total_cloud_cover, flux%lw_derivatives)
    end if

    ! Loop through columns
    do jcol = istartcol,iendcol
      if (total_cloud_cover(jcol) < config%cloud_fraction_threshold) then
        ! No cloud in profile and clear-sky fluxes already
        ! calculated: copy them over
        flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
        flux%lw_dn_surf_g(:,jcol) = flux%lw_dn_surf_clear_g(:,jcol)
      end if
    end do
    if (config%do_lw_derivatives) then
      mask_1d = .FALSE.
      do jcol = istartcol,iendcol
        if (total_cloud_cover(jcol) < config%cloud_fraction_threshold) then
          mask_1d(jcol) = .TRUE.
        end if
      end do
      call calc_lw_derivatives_ica(ng, nlev, istartcol, iendcol, mask_1d, trans_clear(:,:,:), flux_up_clear(:,nlev+1,:), &
           &                       flux%lw_derivatives)
 
    end if ! Cloud is present in profile

    if (lhook) call dr_hook('radiation_mcica_lw_acc:solver_mcica_lw',1,hook_handle)
    
  end subroutine solver_mcica_lw

end module radiation_mcica_lw_acc
