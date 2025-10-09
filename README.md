# Cinemática Inversa del Robot ABB IRB140

Este proyecto implementa y analiza la resolución del problema de cinemática inversa de forma numerica para el robot industrial ABB IRB140, utilizando el modelo cinemático basado en los parámetros Denavit–Hartenberg (DH).

# Descripción general

El trabajo aborda el cálculo de las variables articulares del robot a partir de una pose objetivo (posición y orientación deseadas del efector final), explorando tanto métodos analíticos como numéricos iterativos.
Se busca maximizar la precisión y la tasa de convergencia del algoritmo mediante:

* Generación de semillas iniciales inteligentes basadas en la pose y la configuración deseada.

* Implementación de un método iterativo de control de posición (basado en la cinemática diferencial y la pseudo-inversa amortiguada del Jacobiano).

* Estrategias de envoltura angular (“joint wrapping”) y manejo circular de articulaciones para evitar divergencias.

* Incorporación de un esquema multiintento con perturbaciones para garantizar la convergencia incluso en zonas de singularidad.

# Estructura

Se podrán encontrar los archivos:

* cinematica_irb140.ipynb: que contiene el desarrollo teórico y práctico del proyecto.
* generar_animacion.py: toma la implementacion de la clase creada en el notebook y captura frame por frame la evolucion de la trayectoria del algoritmo. Se generar archivos .png que pueden ser convertidos a video empleando la herramiente ffmpeg
* venv_anim: entorno virtual de Python preconfigurado con las dependencias necesarias para ejecutar los scripts y notebooks.
* frames: carpeta en el que se guardan los archivos .png generados para crear la animacion


Se emplea el siguiente comando para la generacion del video:

```bash
ffmpeg -framerate 60 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p animacion.mp4
```

  