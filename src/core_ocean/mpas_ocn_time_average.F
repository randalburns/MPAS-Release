! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module ocn_time_average

    use mpas_grid_types

    implicit none
    save
    public

    contains 

    subroutine ocn_time_average_init(average)!{{{
        type (average_type), intent(inout) :: average

        real (kind=RKIND), pointer :: nAverage

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: avgVelocityZonal, avgVelocityMeridional, varVelocityZonal, varVelocityMeridional
        real (kind=RKIND), dimension(:,:), pointer :: avgNormalVelocity, varNormalVelocity, avgVertVelocityTop

        real (kind=RKIND), dimension(:,:), pointer :: avgVelocityX, avgTemperature, avgTracer1, avgLayerThickness, avgPlumeVelocityX, avgNPlumeVelocityX

        nAverage => average % nAverage % scalar

        avgSSH                => average % avgSSH % array
        varSSH                => average % varSSH % array
        avgVelocityZonal      => average % avgVelocityZonal % array
        avgVelocityMeridional => average % avgVelocityMeridional % array
        varVelocityZonal      => average % varVelocityZonal % array
        varVelocityMeridional => average % varVelocityMeridional % array
        avgNormalVelocity     => average % avgNormalVelocity % array
        varNormalVelocity     => average % varNormalVelocity % array
        avgVertVelocityTop    => average % avgVertVelocityTop % array

        avgVelocityX    => average % avgVelocityX % array
        avgTemperature    => average % avgTemperature % array
        avgTracer1    => average % avgTracer1 % array
        avgLayerThickness    => average % avgLayerThickness % array
        avgPlumeVelocityX    => average % avgPlumeVelocityX % array
        avgNPlumeVelocityX    => average % avgNPlumeVelocityX % array

        nAverage = 0

        avgSSH = 0.0
        varSSH = 0.0
        avgVelocityZonal = 0.0
        avgVelocityMeridional = 0.0
        varVelocityZonal = 0.0
        varVelocityMeridional = 0.0
        avgNormalVelocity = 0.0
        varNormalVelocity = 0.0
        avgVertVelocityTop = 0.0

        avgVelocityX = 0.0 
        avgTemperature = 0.0 
        avgTracer1 = 0.0 
        avgLayerThickness = 0.0 
        avgPlumeVelocityX = 0.0
        avgNPlumeVelocityX = 0.0

    end subroutine ocn_time_average_init!}}}

    subroutine ocn_time_average_accumulate(average, state, diagnostics)!{{{
        type (average_type), intent(inout) :: average
        type (state_type), intent(in) :: state
        type (diagnostics_type), intent(in) :: diagnostics

        real (kind=RKIND), pointer :: nAverage, old_nAverage

        real (kind=RKIND), dimension(:), pointer :: ssh
        real (kind=RKIND), dimension(:,:), pointer :: normalVelocityZonal, normalVelocityMeridional, normalVelocity, vertVelocityTop

        real (kind=RKIND), dimension(:,:), pointer :: avgNormalVelocity, varNormalVelocity, avgVertVelocityTop
        real (kind=RKIND), dimension(:,:), pointer :: avgVelocityZonal, avgVelocityMeridional, varVelocityZonal, varVelocityMeridional
        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH

        real (kind=RKIND), dimension(:,:), pointer :: old_avgNormalVelocity, old_varNormalVelocity, old_avgVertVelocityTop
        real (kind=RKIND), dimension(:,:), pointer :: old_avgVelocityZonal, old_avgVelocityMeridional, old_varVelocityZonal, old_varVelocityMeridional
        real (kind=RKIND), dimension(:), pointer :: old_avgSSH, old_varSSH

        real (kind=RKIND), dimension(:,:), pointer :: avgVelocityX, avgTemperature, avgTracer1, avgLayerThickness, avgPlumeVelocityX, avgNPlumeVelocityX

        nAverage     => average % nAverage  % scalar

        normalVelocity  => state % normalVelocity % array
        ssh             => state % ssh % array

        normalVelocityZonal      => diagnostics % normalVelocityZonal % array
        normalVelocityMeridional => diagnostics % normalVelocityMeridional % array
        vertVelocityTop          => diagnostics % vertVelocityTop % array

        avgSSH                   => average % avgSSH % array
        varSSH                   => average % varSSH % array
        avgVelocityZonal         => average % avgVelocityZonal % array
        avgVelocityMeridional    => average % avgVelocityMeridional % array
        varVelocityZonal         => average % varVelocityZonal % array
        varVelocityMeridional    => average % varVelocityMeridional % array
        avgNormalVelocity        => average % avgNormalVelocity % array
        varNormalVelocity        => average % varNormalVelocity % array
        avgVertVelocityTop       => average % avgVertVelocityTop % array

        avgVelocityX    => average % avgVelocityX % array
        avgTemperature    => average % avgTemperature % array
        avgTracer1    => average % avgTracer1 % array
        avgLayerThickness    => average % avgLayerThickness % array
        avgPlumeVelocityX    => average % avgPlumeVelocityX % array
        avgNPlumeVelocityX    => average % avgNPlumeVelocityX % array

        avgSSH = avgSSH + ssh
        varSSH = varSSH + ssh**2
        avgVelocityZonal = avgVelocityZonal + normalVelocityZonal
        avgVelocityMeridional = avgVelocityMeridional + normalVelocityMeridional
        varVelocityZonal = varVelocityZonal + normalVelocityZonal**2
        varVelocityMeridional = varVelocityMeridional + normalVelocityMeridional**2
        avgNormalVelocity = avgNormalVelocity + normalVelocity
        varNormalVelocity = varNormalVelocity + normalVelocity**2
        avgVertVelocityTop = avgVertVelocityTop + vertVelocityTop

        avgVelocityX    = avgVelocityX + diagnostics % normalVelocityX % array
        avgTemperature    = avgTemperature + state % tracers % array(1,:,:)
        avgTracer1    = avgTracer1 + state % tracers % array(3,:,:)
        avgLayerThickness    = avgLayerThickness + state % layerThickness % array

        where(state % tracers % array(3,:,:).gt.0.01)
           avgPlumeVelocityX    = avgPlumeVelocityX + diagnostics % normalVelocityX % array
           avgNPlumeVelocityX    = avgNPlumeVelocityX + 1.0
        end where

        nAverage = nAverage + 1
    end subroutine ocn_time_average_accumulate!}}}

    subroutine ocn_time_average_normalize(average)!{{{
        type (average_type), intent(inout) :: average

        real (kind=RKIND), pointer :: nAverage

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: avgVelocityZonal, avgVelocityMeridional, varVelocityZonal, varVelocityMeridional
        real (kind=RKIND), dimension(:,:), pointer :: avgNormalVelocity, varNormalVelocity, avgVertVelocityTop

        real (kind=RKIND), dimension(:,:), pointer :: avgVelocityX, avgTemperature, avgTracer1, avgLayerThickness, avgPlumeVelocityX, avgNPlumeVelocityX

        nAverage => average % nAverage  % scalar

        avgSSH                => average % avgSSH % array
        varSSH                => average % varSSH % array
        avgVelocityZonal      => average % avgVelocityZonal % array
        avgVelocityMeridional => average % avgVelocityMeridional % array
        varVelocityZonal      => average % varVelocityZonal % array
        varVelocityMeridional => average % varVelocityMeridional % array
        avgNormalVelocity     => average % avgNormalVelocity % array
        varNormalVelocity     => average % varNormalVelocity % array
        avgVertVelocityTop    => average % avgVertVelocityTop % array

        avgVelocityX    => average % avgVelocityX % array
        avgTemperature    => average % avgTemperature % array
        avgTracer1    => average % avgTracer1 % array
        avgLayerThickness    => average % avgLayerThickness % array
        avgPlumeVelocityX    => average % avgPlumeVelocityX % array
        avgNPlumeVelocityX    => average % avgNPlumeVelocityX % array

        if(nAverage > 0) then
          avgSSH = avgSSH / nAverage
          varSSH = varSSH / nAverage
          avgVelocityZonal = avgVelocityZonal / nAverage
          avgVelocityMeridional = avgVelocityMeridional / nAverage
          varVelocityZonal = varVelocityZonal / nAverage
          varVelocityMeridional = varVelocityMeridional / nAverage
          avgNormalVelocity = avgNormalVelocity / nAverage
          varNormalVelocity = varNormalVelocity / nAverage
          avgVertVelocityTop = avgVertVelocityTop / nAverage

          avgVelocityX    = avgVelocityX / nAverage
          avgTemperature    = avgTemperature / nAverage
          avgTracer1    = avgTracer1 / nAverage
          avgLayerThickness    = avgLayerThickness / nAverage
          avgPlumeVelocityX    = avgPlumeVelocityX / avgNPlumeVelocityX

          nAverage = 0
        end if
    end subroutine ocn_time_average_normalize!}}}

end module ocn_time_average
