! A fortran95 program for G95
! By GRUPO 3 Informatica 1ºGIA ETSIAE: Manuel, Vito, Victor, David y Adrian
program main

  use utilidades

  implicit none

  !_________________________________________________________________________________________________
  !DECLARACION E INICIALIZACION DE VARIABLES DE PARTIDA

  !Declaracion de variables de partida
  integer, parameter:: nPlazas_max=180, nClases=5
  real*8:: media(1:nClases), precio(1:nClases), desvest(1:nClases)


  !Declaracion de variables intermedias, previas al resultado final
  real*8:: emsr(nClases,nPlazas_max), maxEmsr(2,nPlazas_max)


  !Declaracion de variables finales
  integer:: nPlazasFin(nClases)
  real*8:: beneficios


  !Inicializacion de variables de partida
  media(2:nClases) = reshape((/9d0,23d0,44d0,66d0/),(/nClases-1/))
  !media(1) => se resuelve por sist. ecuaciones


  desvest(2:nClases) = reshape((/4d0,7d0,9d0,15d0/), (/nClases-1/))
  !desvest(1) => se resuelve por sist. ecuaciones


  precio = reshape((/180d0,140d0,100d0,70d0,35d0/), (/nClases/))

  !____________________________________________ENCABEZADO_______________________________________________

  print*, "---------------------------------------------------------------------------------"
  print*, "                               GRUPO 3 Informatica                               "
  print*, "---------------------------------------------------------------------------------"
  print*, ""

  !_____________________________________________________________________________________________________
  !1º) RESOLUCION CLASE A

  print*,      " [PRODUCTO A]"
  print*,      ""
  print*,      " ______________________________________","     ---------"
  call newton(media,desvest,nClases)
  write(*,100) "*La media de la clase A es:            ", media(1)
  write(*,100) "*La desviacion tipica de la clase A es:", desvest(1)
  100 format(A40,2x,F11.6)
  print*,      " ______________________________________","     ---------"
  print*,      ""
  print*,      ""


  !______________________________________________________________________________________________________
  !2º) CALCULO EMSR POR TRAPECIO
  call trapecio(emsr,desvest,media,precio,nClases,nPlazas_max)


  !______________________________________________________________________________________________________
  !3º) ORDENACION DE LOS EMSR
  maxEmsr = max180(emsr,nClases,nPlazas_max)

  !Datos de chequeo
  print*, " [VALORES EMSR ORDENADOS]"
  print*, ""
  print*, "___________________         -----------------------------------------------------"
  write(*,200) "EMSR NUMERO 7 ", maxEmsr(1,7)
  print*, "___________________         -----------------------------------------------------"
  write(*,200) "EMSR NUMERO 19", maxEmsr(1,19)
  200 format(A15,30x,F11.6)
  print*, "___________________         -----------------------------------------------------"
  print*, ""
  print*, ""


  !_______________________________________________________________________________________________________
  !4º) RECUENTO DE LOS EMSR
  nPlazasFin = recuento(int(maxEmsr(2,:)),nPlazas_max,nClases)   !Recuento distribucion final de asientos


  !_______________________________________________________________________________________________________
  !5º) RESULTADOS OBTENIDOS
  print*," [DISTRIBUCION DE ASIENTOS]"
  print*,""
  print*,"___________________         -----------------------------------------------------"
  print*,"PRODUCTO                     /A/       /B/         /C/         /D/         /E/ "
  print*,"___________________         -----------------------------------------------------"
  print*,"NUMERO DE ASIENTOS",nPlazasFin
  print*,"___________________         -----------------------------------------------------"
  beneficios = income (nPlazasFin,precio,nClases)
  write(*,300) "INGRESOS OPTIMIZADOS", beneficios, "$"
  300 format(A21,25x,F12.2,A2)
  print*,"___________________         -----------------------------------------------------"
  print*,""
  print*,"______________________________________________"
  print*,"Autores: Manuel , Vito, Victor, David y Adrian"


end
