# Proyecto EDP Fortran

Este proyecto implementa la resolución numérica de una ecuación en derivadas parciales (EDP) usando el método de Crank-Nicolson en Fortran. Incluye dos experimentos principales para analizar el comportamiento del método bajo diferentes condiciones.

## Estructura del proyecto

- **experiment_A.f90**: Ejecuta el Experimento A, donde el paso temporal `k` es fijo.
- **experiment_B.f90**: Ejecuta el Experimento B, donde el parámetro `λ` es fijo.
- **utils_solver.f90**: Contiene el módulo `utils_solver` con el método de Thomas para resolver sistemas tridiagonales.
- **.gitignore**: Archivos y ejecutables ignorados por git.

## Compilación

Desde la terminal, ejecuta:

```bash
# Compilar tus experimentos originales
gfortran -c utils_solver.f90
gfortran utils_solver.o experiment_A.f90 -o experiment_A
gfortran utils_solver.o experiment_B.f90 -o experiment_B

# Compilar nueva implementación del método explícito 
gfortran -O3 -o heat_fortran calor_1d.f90
g++ -std=c++17 -O3 -o heat_cpp calor_1d.cc

# Ejecutar
./experiment_A
./experiment_B
./heat_fortran
./heat_cpp
```


## Notas

- Asegúrate de tener instalado `gfortran`.
- Los archivos `.o` y `.mod` generados por la compilación están ignorados en el repositorio.
- Puedes modificar los parámetros en los archivos de los experimentos para probar diferentes configuraciones.

---
Proyecto para prácticas de métodos numéricos en Fortran.