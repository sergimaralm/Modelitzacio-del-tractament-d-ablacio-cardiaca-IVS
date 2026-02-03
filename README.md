# Modelización del Tratamiento de Ablación Cardíaca (IVS) 🫀🔥

Este repositorio contiene la implementación de diversos métodos numéricos para simular la evolución de la temperatura en el Tabique Interventricular (IVS) durante un tratamiento de ablación cardíaca.

El objetivo principal es resolver la **Ecuación del Calor** unidimensional para modelar cómo se propaga la energía térmica en el tejido, comparando diferentes esquemas de discretización temporal.

## 📋 Contenido del Repositorio

El proyecto está estructurado en módulos de Python, cada uno implementando un algoritmo numérico específico:

| Archivo | Descripción |
| :--- | :--- |
| `eulerexplicito.py` | Implementación del método de **Euler Explícito**. |
| `eulerimplicit.py` | Implementación del método de **Euler Implícito** (incondicionalmente estable). |
| `crank-nicholson.py` | Implementación del método de **Crank-Nicolson** (orden 2, incondicionalmente estable). |
| `soluciotemps.py` | Análisis de la evolución temporal y comparativas. |
| `funcions.py` | Librería auxiliar con funciones compartidas (solución analítica, método de Jacobi, etc.). |
| `figures/` | Carpeta donde se guardan las gráficas generadas automáticamente. |
| `Pràctica 1.pdf` | Documentación teórica y enunciado de la práctica. |

## 🚀 Requisitos e Instalación

Para ejecutar este proyecto, necesitarás **Python 3.x** y las librerías estándar de cálculo científico.

1. Clona el repositorio:
   ```bash
   git clone [https://github.com/sergimaralm/Modelitzacio-del-tractament-d-ablacio-cardiaca-IVS.git](https://github.com/sergimaralm/Modelitzacio-del-tractament-d-ablacio-cardiaca-IVS.git)
   cd Modelitzacio-del-tractament-d-ablacio-cardiaca-IVS
