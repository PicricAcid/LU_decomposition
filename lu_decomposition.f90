program LU_decomposition
  implicit none
!Parameters
  integer, parameter:: n = 4
!variables
  integer:: i, j, k, unit, imax, jmax, t
  double precision:: sum, cir1, cir3, pivot
!arrays
  double precision, allocatable :: a(:, :), l(:, :), u(:, :), b(:), x(:), y(:), jp(:)
!init variables
  allocate(a(n, n), l(n, n), u(n, n), b(n), x(n), y(n), jp(n))
  l = 0.0d0
  u = 0.0d0

  unit = 10
  open(unit, file='array.csv', status='old')
  read(unit, '()')
  do i=1, n
    read(unit, *) (a(i, j), j=1, n)
  end do
  read(unit, '()')
  read(unit, *) (b(i), i=1, n)
  close(unit)

  do i=1, n
    jp(i) = i
  end do

!solve the equations Ax = b for x
  do k=1, n
!pivotting
    pivot = dabs(a(k, k))
    imax = k
    jmax = k
    do i=k, n
      do j=k, n
        if(dabs(a(i, k)) > pivot) then
          pivot = dabs(a(i, k))
          imax = i
          jmax = j
        end if
      end do
    end do

    if(jmax /= k) then
      t = jp(k)
      jp(k) = jmax
      jp(jmax) = t
    end if

    do i=1, n
      pivot = a(i, k)
      a(i, k) = a(i, jmax)
      a(i, jmax) = pivot
    end do

    do j=1, n
      pivot = a(k, j)
      a(k, j) = a(imax, j)
      a(imax, j) = pivot
    end do

    pivot = b(k)
    b(k) = b(imax)
    b(imax) = pivot

!LU decomposition
    cir1 = 1.0d0/a(k, k)
    u(k, k) = a(k, k) !cir1
    do i=k+1, n
      a(i, k) = a(i, k)*cir1
    end do
    do j=k+1, n
      cir3 = a(k, j)
      u(k, j) = cir3
      do i=k+1, n
        a(i, j) = a(i, j) - a(i, k)*cir3
      end do
    end do
  end do

  do i=1, n
    do j=1, n
      if(i == j) then
        l(i, j) = 1.0d0
      else if(i > j) then
        l(i, j) = a(i, j)
      end if
    end do
  end do

! forward substitution
  do i=1, n
    sum = 0
    do k=1, i-1
      sum = sum + l(i, k)*y(k)
    end do
    y(i) = b(i) - sum
  end do

!back substitution
  do i=n, 1, -1
    sum = 0
    do k=i+1, n
      sum = sum + u(i, k)*x(k)
    end do
    x(i) = (y(i) - sum)/u(i, i)
  end do

  do i=1, n
    write(*, *) (a(i, j), j=1, n), b(i)
  end do

  do i=1, n
    write(*, *) x(i)
  end do

end program
