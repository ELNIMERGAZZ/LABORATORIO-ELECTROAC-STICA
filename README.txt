# Simulador mecánico 5 GDL — v4

## Topología implementada

El modelo usa cinco grados de libertad:

q = [u₁, u₂, u₃, u₄, u₅]ᵀ

La topología sigue el esquema mecánico:

Pared — F(t) — M₁ — K₁ — nodo sin masa

Desde ese nodo salen dos ramas paralelas:

- Rama superior: K₂ — M₂ — B₂
- Rama inferior: K₃ — M₃ — B₃

Las dos ramas se reúnen y llegan a M₄.

Después:

- B₄ entre M₄ y M₅
- B₅ entre M₅ y el piso
- K₄ entre M₅ y la pared derecha
- B₁ entre M₁ y el piso

El nodo entre K₁, K₂ y K₃ no tiene masa.

## Cambios principales

- Slider y campo numérico sincronizados en ambos sentidos.
- Los sliders actualizan inmediatamente el modelo.
- La simulación tiene una velocidad predeterminada menor: 0.35×.
- Las velocidades se muestran con formato decimal legible, sin notación e-.
- El panel de velocidades está debajo del diagrama.
- El diagrama ocupa mucho más espacio horizontal y vertical.
- Se eliminó la parte de espectro, movilidad y coherencia.
- Se conservó la respuesta temporal, la telemetría y la exportación CSV.
- Se añadió una velocidad promedio calculada a partir de las velocidades absolutas registradas de las cinco masas.
- Se retiraron de la interfaz los paneles de dinámica modal y balance energético.

## Ejecución

Abrir `index.html` directamente en Chrome, Edge o Firefox.
No requiere servidor, Python ni librerías externas.
