module spline_kernal

  implicit none

contains

  subroutine w(q, w_val, w_grad)
    real, intent(in) :: q
    real :: w_val, w_grad, sigma

    !normalisation
    sigma = (2./3.)


    !peicewise implementation of the spline kernal from Price, D., 2010
    if (q < 1) then
      w_val = (1./4.)*(2-q)**(3) - (1-q)**(3)
      w_grad = 3.*(1-q)**(2) - (3./4.)*(2-q)**(2)
    else if (q < 2) then
      w_val = (1./4.)*(2-q)**(3)
      w_grad = -(3./4.)*(2-q)**(2)
    else
      w_val = 0
      w_grad = 0
    end if

    w_val = sigma * w_val
    w_grad = sigma * w_grad

  end subroutine W



end module spline_kernal
