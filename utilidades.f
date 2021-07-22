       module utilidades
            contains
!___________________________________________________________________________
!1º) RESOLUCION CLASE A

            !Resolución de desvest(1) y media(1) en la clase A
            subroutine newton (c,d,m)
                integer, intent(in):: m
                real*8, intent(inout):: c(m), d(m)

                integer,parameter:: p = 2, nmax = 10
                real*8,parameter:: tol = 1d-6                                 !Error de la solución ajustable

                integer:: i, n
                real*8:: x(p), y(p), a(p,p), errf, errx

                x(1) = 1.6112293                                              !Valores iniciales: primera aproximacion, dados por Taylor de la ecuacion logaritmica centrado en 2
                x(2) = 5.0151694                                              !(VER DOCUMENTO GRAFICO). UTILIZO 7 DECIMALES
                n = 0

                !Terminos independientes para el primer sistema
                call terminos_independientes(x,y,p)                           !Los term. indep. son el valor negativo de las funciones dadas igualadas a 0 (VER METODO DE NEWTON)

                do n= 1, nmax

                    call jacobiano(x,a,p)                                     !Metemos jacobiano en la matriz a
                                                                              !Aplicamos LU al sistema: jacobiano*dx=terminos independientes
                    call factorizar(a,p)
                    call sustituir(a,y,p)                                     !Una vez aqui, tenemos en "y" la solucion del sistema anterior, es decir, el valor de "dx"
                                                                              !Ahora tenemos en cuenta que la segunda aproximacion (x2) serà de la forma:
                                                                              !x2=x1+dx
                    do i= 1, p
                        x(i) = x(i)+y(i)
                    end do


                    errx = norma(y,p)/norma(x,p)                              !Calculo de los dos errores a tener en cuenta(errx y errf)--> ver metodo de Newton
                    call terminos_independientes(x,y,p)
                    errf = norma(y,p)
                    if (max(errx,errf) .lt. tol) exit
                end do
                c(1) = x(1)
                d(1) = x(2)

            end subroutine newton


            !Term. indp. sistema de ecuaciones no lineales
            subroutine terminos_independientes( x,y,p)
                implicit none
                integer, intent(in):: p
                real*8,intent(in):: x(p)
                real*8,intent(inout):: y(p)

                    y(1) = -exp(x(1))-log(x(1))+x(2)
                    y(2) = (x(1))**2d0-(x(1))-6d0+x(2)

            end subroutine terminos_independientes


            !Jacobiano sistema de ecuaciones no lineales
            subroutine jacobiano(x,a,p)
                implicit none
                integer, intent(in):: p
                real*8, intent(in):: x(p)
                real*8, intent(inout):: a(p,p)

                    a(1,1) = exp(x(1))+1d0/(x(1))
                    a(1,2) = -1d0
                    a(2,1) = -2d0*(x(1))+1d0
                    a(2,2) = -1d0

            end subroutine jacobiano


            !Funcion norma euclidea
            function norma (a,n)
                implicit none
                integer:: n, i
                real*8:: a(n), norma
                norma = 0d0

                do i= 1, n
                    norma = norma+a(i)**2d0
                end do
                norma = sqrt(norma)

            end function norma


            !Factorizacion por el metodo LU
            subroutine factorizar (a,n)
                implicit none
                integer, intent(in):: n
                real*8,intent(inout):: a(n,n)

                integer :: k, i, h, j
                real*8 :: s

                do k = 1, n
                    do i = k, n
                    s = 0d0
                        do h=1, k-1
                        s = s + a(i,h) * a(h,k)
                        end do
                        a(i,k) = a(i,k) - s
                    end do
                    do j = k+1, n
                        s = 0d0
                        do h = 1, k-1
                            s = s + a(k,h) * a(h,j)
                        end do
                    a(k,j) = (1d0 / a(k,k)) * (a(k,j) - s)
                    end do
                end do

            end subroutine factorizar


            !Sustuir valores de matrices superior e inferior factorizadas por LU
            subroutine sustituir (a,b,n)
            implicit none
                integer, intent(in) :: n
                real*8, intent(inout) :: a(n,n), b(n)

                integer :: k, h
                real*8:: s

                    do k = 1, n
                        s = 0d0
                        do h = 1, k-1
                            s = s + a(k,h) * b(h)
                        end do
                            b(k) = (b(k) - s) / a(k,k)
                    end do
                do k = n, 1, -1
                    s = 0d0
                    do h = k+1, n
                        s = s+a(k,h)*b(h)
                    end do
                    b(k) = b(k) - s
                end do

            end subroutine sustituir


!___________________________________________________________________________________________
!2º) CALCULO EMSR POR TRAPECIO

            !Funcion de una distribucion de probabilidad normal, uso en trapecio
            function normal (x,media,desvest)
                implicit none
                real*8::x, normal, desvest, media, c

                c = (1d0/(desvest*sqrt(2d0*acos(-1d0))))
                normal = c*exp(-((x-media)**2d0)/(2d0*(desvest**2d0)))

            end function normal


           !Integral para calcular los emsr
          subroutine trapecio(emsr,desvest,media,precio,nClases,nPlazas)
                implicit none
                real*8,intent(in) :: desvest(nClases), media(nClases)
                real*8,intent(in) :: precio(nClases)
                integer, intent(in) :: nClases, nPlazas
                real*8,intent(inout) :: emsr(nClases,nPlazas)

                !Parametros ajustables para los limites e intervalos de la integral
                integer, parameter:: iter_trapecio=500
                real*8,parameter:: limit_sup=300d0

                integer :: i, k, j
                real*8 :: x, h, t, s, m



                do k= 1, nClases                                            !Itera para cada una de las cinco clases
                    do j= 1, nPlazas                                        !Itera para las 180 primeras plazas de cada clase
                        h = (limit_sup-dble(j))/dble(iter_trapecio)         !Amplitud h del intervalo a integrar
                        s = 0d0
                        do i= 1, iter_trapecio-1
                            x = dble(j)+h*dble(i)                           !Se recorre todo el dominio de integracion [j,limit_sup] sumando los valores de la funcion
                            s = s+normal(x,media(k),desvest(k))
                        end do

                       m = (0.5*normal(limit_sup,media(k),desvest(k))+s)
                       t = h*(0.5*normal(dble(j),media(k),desvest(k))+m)    !Aplicacion formula trapecio

                    emsr(k,j) = t*precio(k)                                 !EMSR=probabilidad*precio, almacenados en matriz emsr(nClases,nPlazas)
                    enddo
                enddo

            end subroutine trapecio


!_________________________________________________________________________
!3º) ORDENACION DE LOS EMSR

            !Calcula los 180 mayores esmr y los asocia con su clase (Clase A =  1, ... , clase E = 5)
            function max180 (emsr,nClases,nPlazas_max)
                real*8,intent(in):: emsr(nClases,nPlazas_max)
                integer,intent(in):: nClases, nPlazas_max

                real*8 :: aux(nClases,nPlazas_max)
                real*8:: max180(2,nPlazas_max), mVal                        !max180: matriz que contiene en la primera fila los 180 mayores ESMR y en la segunda su clase asociada
                integer :: i, j, k, r, s

                aux = emsr

                do k = 1, nPlazas_max
                    mVal = 0d0
                    do i = 1, nClases                                       !Recorremos toda la matriz emsr(nClases,nPlazas_max)
                        do j = 1, nPlazas_max
                            if (aux(i,j) > mVal) then                       !Comparamos cada valor de la matriz emsr(nClases,nPlazas_max) con el resto

                                mVal = aux(i,j)                             !Guardamos el valor del mayor valor encontrado
                                r = i                                       !Guardamos en r,s la posición de dicho mayor valor encontrado
                                s = j

                            end if
                        end do
                    end do
                    max180(1,k) = mVal                                      !Colocamos en la fila 1 el mayor valor ordenado
                    max180(2,k) = dble(r)                                   !Colocamos en la fila 2 la clase asociada con el valor encontrado
                    aux(r,s) = 0d0                                          !El mayor valor encontrado lo igualamos a 0 para no volver a seleccionarlo en sucesivos bucles
                end do

            end function max180


!_________________________________________________________________________
!4º) RECUENTO DE LOS EMSR

            function recuento (v,num_asientos,num_clases)
                implicit none
                integer, intent(in) :: num_asientos, v(num_asientos)         !Vector que contiene las clases de cada asiento: Clase A = 1,..., clase E = 5.
                integer, intent(in) :: num_clases

                integer :: recuento(num_clases), i

                recuento = 0                                                 !Comenzamos el recuento en 0
                do i = 1, num_asientos
                    recuento(v(i)) = recuento(v(i)) + 1                      !Se suma uno en el vector recuento(num_clases) en la posicion correspondiente a la clase que sea
                end do                                                       !Recuento(1) contendra el nº de la clase A, recuento(2) el de la clase(2)...

            end function recuento


!_________________________________________________________________________
!5º) RESULTADOS OBTENIDOS

            function income (nPlazasFin,precio,nClases)                      !Funcion que calcula los beneficios obtenidos con la distribucion optima de plazas
                implicit none
                integer, intent(in):: nPlazasFin(nClases), nClases
                real*8,intent(in):: precio(nClases)

                real*8:: income
                integer::i

                income = 0d0

                do i= 1, nClases
                    income = income+precio(i)*nPlazasFin(i)
                end do

            end function



        end module utilidades
