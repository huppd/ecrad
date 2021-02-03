! radiation_adding_ica_lw_acc.F90 - Longwave adding method in independent column approximation
!
! Copyright (C) 2015-2017 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-07-12  R. Hogan  Fast adding method for if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_adding_ica_lw_acc

  public

contains

  !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering, by successively adding the contribution of
  ! layers starting from the surface to compute the total albedo and
  ! total upward emission of the increasingly larger block of
  ! atmospheric layers.
  subroutine adding_ica_lw(ng, nlev, istartcol, iendcol, mask, &
       &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
       &  flux_up, flux_dn)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs
    integer, intent(in) :: ng ! number of spectral intervals
    integer, intent(in) :: nlev, istartcol, iendcol ! number of levels

    ! where to do this computation
    logical, intent(in) :: mask(istartcol:iendcol)

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng, nlev)   :: reflectance, transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng, nlev)   :: source_up, source_dn

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(istartcol:iendcol, ng, nlev+1) :: flux_up, flux_dn
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), dimension(istartcol:iendcol, ng, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to emission below
    ! that half-level (W m-2)
    real(jprb), dimension(istartcol:iendcol, ng, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), dimension(istartcol:iendcol, ng, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jg, jlev, jcol

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_lw_acc:adding_ica_lw',0,hook_handle)

    albedo(:,nlev+1,:) = albedo_surf

    ! At the surface, the source is thermal emission
    source(:,nlev+1,:) = emission_surf

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    ! Loop through columns
    do jlev = nlev,1,-1
      do jcol = istartcol,iendcol
      ! Next loop over columns. We could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  Rather, we do it with an explicit
      ! loop.
        do jg = 1,ng
          if (mask(jcol)) then
            ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
            inv_denominator(jcol,jg,jlev) = 1.0_jprb &
                 &  / (1.0_jprb-albedo(jcol,jg,jlev+1)*reflectance(jcol,jg,jlev))
            ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
            albedo(jcol,jg,jlev) = reflectance(jcol,jg,jlev) + transmittance(jcol,jg,jlev)*transmittance(jcol,jg,jlev) &
                 &  * albedo(jcol,jg,jlev+1) * inv_denominator(jcol,jg,jlev)
            ! Shonk & Hogan (2008) Eq 11:
            source(jcol,jg,jlev) = source_up(jcol,jg,jlev) &
                 &  + transmittance(jcol,jg,jlev) * (source(jcol,jg,jlev+1) &
                 &                    + albedo(jcol,jg,jlev+1)*source_dn(jcol,jg,jlev)) &
                 &                   * inv_denominator(jcol,jg,jlev)
          end if
        end do
      end do
    end do

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(:,1,:) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to emission
    ! below that level
    flux_up(:,1,:) = source(:,1,:)

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
    do jlev = 1,nlev
      ! Loop through columns
      do jcol = istartcol,iendcol
        do jg = 1,ng
          if (mask(jcol)) then
            ! Shonk & Hogan (2008) Eq 14 (after simplification):
            flux_dn(jcol,jg,jlev+1) &
                 &  = (transmittance(jcol,jg,jlev)*flux_dn(jcol,jg,jlev) &
                 &     + reflectance(jcol,jg,jlev)*source(jcol,jg,jlev+1) &
                 &     + source_dn(jcol,jg,jlev)) * inv_denominator(jcol,jg,jlev)
            ! Shonk & Hogan (2008) Eq 12:
            flux_up(jcol,jg,jlev+1) = albedo(jcol,jg,jlev+1)*flux_dn(jcol,jg,jlev+1) &
                 &            + source(jcol,jg,jlev+1)
          end if
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_adding_ica_lw_acc:adding_ica_lw',1,hook_handle)

  end subroutine adding_ica_lw


  !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering in cloudy layers only.
  subroutine fast_adding_ica_lw(ng, nlev, istartcol, iendcol, mask, &
       &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
       &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
       &  flux_up, flux_dn)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs
    integer, intent(in) :: ng ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels
    integer, intent(in) :: istartcol
    integer, intent(in) :: iendcol

    logical, intent(in) :: mask(istartcol:iendcol)

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng, nlev)   :: reflectance, transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng, nlev)   :: source_up, source_dn

    ! Determine which layers are cloud-free
    logical, intent(in) :: is_clear_sky_layer(istartcol:iendcol, nlev)

    ! Index to highest cloudy layer
    integer, intent(in) :: i_cloud_top(istartcol:iendcol)

    ! Pre-computed clear-sky downwelling fluxes (W m-2) at half-levels
    real(jprb), intent(in), dimension(istartcol:iendcol, ng, nlev+1)  :: flux_dn_clear

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(istartcol:iendcol, ng, nlev+1) :: flux_up, flux_dn
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), dimension(istartcol:iendcol, ng, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to emission below
    ! that half-level (W m-2)
    real(jprb), dimension(istartcol:iendcol, ng, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), dimension(istartcol:iendcol, ng, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jlev, jcol, jg

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_lw_acc:fast_adding_ica_lw',0,hook_handle)

    do jcol=istartcol,iendcol
      if (mask(jcol)) then
        ! Copy over downwelling fluxes above cloud from clear sky
        flux_dn(jcol,:,1:i_cloud_top(jcol)) = flux_dn_clear(jcol,:,1:i_cloud_top(jcol))

        albedo(jcol,:,nlev+1) = albedo_surf(jcol,:)

        ! At the surface, the source is thermal emission
        source(jcol,:,nlev+1) = emission_surf(jcol,:)
      endif
    end do

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    ! Reflectance of this layer is zero, simplifying the expression
    do jcol = istartcol,iendcol
      do jg = 1,ng
        do jlev = nlev,i_cloud_top(jcol),-1
          if (mask(jcol)) then
            if (is_clear_sky_layer(jcol,jlev)) then
              albedo(jcol,jg,jlev) = transmittance(jcol,jg,jlev)*transmittance(jcol,jg,jlev)*albedo(jcol,jg,jlev+1)
              source(jcol,jg,jlev) = source_up(jcol,jg,jlev) &
                   &  + transmittance(jcol,jg,jlev) * (source(jcol,jg,jlev+1) &
                   &                    + albedo(jcol,jg,jlev+1)*source_dn(jcol,jg,jlev))
            else
              ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
              inv_denominator(jcol,jg,jlev) = 1.0_jprb &
                   &  / (1.0_jprb-albedo(jcol,jg,jlev+1)*reflectance(jcol,jg,jlev))
              ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
              albedo(jcol,jg,jlev) = reflectance(jcol,jg,jlev) + transmittance(jcol,jg,jlev)*transmittance(jcol,jg,jlev) &
                   &  * albedo(jcol,jg,jlev+1) * inv_denominator(jcol,jg,jlev)
              ! Shonk & Hogan (2008) Eq 11:
              source(jcol,jg,jlev) = source_up(jcol,jg,jlev) &
                   &  + transmittance(jcol,jg,jlev) * (source(jcol,jg,jlev+1) &
                   &                    + albedo(jcol,jg,jlev+1)*source_dn(jcol,jg,jlev)) &
                   &                   * inv_denominator(jcol,jg,jlev)
            end if
          end if
        end do
      end do
    end do

    ! Compute the fluxes above the highest cloud

    do jcol=istartcol,iendcol
      if (mask(jcol)) then
        flux_up(jcol,:,i_cloud_top(jcol)) = source(jcol,:,i_cloud_top(jcol)) &
             &                 + albedo(jcol,:,i_cloud_top(jcol))*flux_dn(jcol,:,i_cloud_top(jcol))
      end if
    end do
    do jcol=istartcol,iendcol
      do jlev = i_cloud_top(jcol)-1,1,-1
        if (mask(jcol)) then
          flux_up(jcol,:,jlev) = transmittance(jcol,:,jlev)*flux_up(jcol,:,jlev+1) + source_up(jcol,:,jlev)
        end if
      end do
    end do

    ! Work back down through the atmosphere from cloud top computing
    ! the fluxes at each half-level
    do jcol=istartcol,iendcol
      do jg = 1,ng
        do jlev = i_cloud_top(jcol),nlev
          if (mask(jcol)) then
            if (is_clear_sky_layer(jcol,jlev)) then
              flux_dn(jcol,jg,jlev+1) = transmittance(jcol,jg,jlev)*flux_dn(jcol,jg,jlev) &
                   &               + source_dn(jcol,jg,jlev)
              flux_up(jcol,jg,jlev+1) = albedo(jcol,jg,jlev+1)*flux_dn(jcol,jg,jlev+1) &
                   &               + source(jcol,jg,jlev+1)
            else
              ! Shonk & Hogan (2008) Eq 14 (after simplification):
              flux_dn(jcol,jg,jlev+1) &
                   &  = (transmittance(jcol,jg,jlev)*flux_dn(jcol,jg,jlev) &
                   &     + reflectance(jcol,jg,jlev)*source(jcol,jg,jlev+1) &
                   &     + source_dn(jcol,jg,jlev)) * inv_denominator(jcol,jg,jlev)
              ! Shonk & Hogan (2008) Eq 12:
              flux_up(jcol,jg,jlev+1) = albedo(jcol,jg,jlev+1)*flux_dn(jcol,jg,jlev+1) &
                   &               + source(jcol,jg,jlev+1)
            end if
          end if
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_adding_ica_lw_acc:fast_adding_ica_lw',1,hook_handle)

  end subroutine fast_adding_ica_lw


  !---------------------------------------------------------------------
  ! If there is no scattering then fluxes may be computed simply by
  ! passing down through the atmosphere computing the downwelling
  ! fluxes from the transmission and emission of each layer, and then
  ! passing back up through the atmosphere to compute the upwelling
  ! fluxes in the same way.
  subroutine calc_fluxes_no_scattering_lw(ng, nlev, istartcol, iendcol, mask, &
       &  transmittance, source_up, source_dn, emission_surf, albedo_surf, flux_up, flux_dn)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs
    integer, intent(in) :: ng ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels
    integer, intent(in) :: istartcol
    integer, intent(in) :: iendcol

    logical, intent(in) :: mask(istartcol:iendcol)

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng, nlev)   :: transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng, nlev)   :: source_up, source_dn

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(istartcol:iendcol, ng, nlev+1) :: flux_up, flux_dn
    
    ! Loop index for model level
    integer :: jlev, jg, jcol

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_lw_acc:calc_fluxes_no_scattering_lw',0,hook_handle)

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(:,:,1) = 0.0_jprb

    ! Work down through the atmosphere computing the downward fluxes
    ! at each half-level
!NEC$ outerloop_unroll(8)
    do jlev = 1,nlev
      do jcol=istartcol,iendcol
        do jg = 1,ng
          if (mask(jcol)) then
            flux_dn(jcol,jg,jlev+1) = transmittance(jcol,jg,jlev)*flux_dn(jcol,jg,jlev) + source_dn(jcol,jg,jlev)
          end if
        end do
      end do
    end do

    ! Surface reflection and emission
    flux_up(:,:,nlev+1) = emission_surf(:,:) + albedo_surf(:,:) * flux_dn(:,:,nlev+1)

    ! Work back up through the atmosphere computing the upward fluxes
    ! at each half-level
!NEC$ outerloop_unroll(8)
    do jlev = nlev,1,-1
      do jcol=istartcol,iendcol
        do jg = 1,ng
          if (mask(jcol)) then
            flux_up(jcol,jg,jlev) = transmittance(jcol,jg,jlev)*flux_up(jcol,jg,jlev+1) + source_up(jcol,jg,jlev)
          end if
        end do
      end do
    end do
    
    if (lhook) call dr_hook('radiation_adding_ica_lw_acc:calc_fluxes_no_scattering_lw',1,hook_handle)

  end subroutine calc_fluxes_no_scattering_lw

end module radiation_adding_ica_lw_acc
