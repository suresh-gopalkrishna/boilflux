# Materials

Properties are piecewise-linear tables in temperature, held constant outside the table range.

## Built-in tables

| Name | Range [°C] | Notes |
|---|---|---|
| `AA6082` | 20–600 | aluminium alloy 6082; \(k\) ≈ 170–190 W/(m K) |
| `Nickel` | 20–900 | pure nickel; the \(c_p\) peak and \(k\) minimum near the Curie point (358 °C) are included |
| `C45` | 20–900 | plain carbon steel; ferrite–austenite effects are smeared into \(c_p\) |

!!! warning "Use your own properties for experiments"
    The built-in tables are representative handbook values meant for synthetic studies and
    examples. Property errors pass directly into the recovered flux, roughly in proportion to the
    error in \(\rho c_p\) for the transient part. For quantitative work, measure the properties of
    your own specimens or take them from a trusted source and load them as shown below.

## Custom tables

```python
mat = bf.Material.from_csv("my_alloy.csv")     # columns: T_C,k,rho,cp (with header)
mat = bf.Material(name="X", T_table=[20, 400, 800], k_table=[45, 38, 27],
                  rho_table=[7850, 7720, 7600], cp_table=[460, 600, 650])
mat = bf.Material.constant(k=50.0, rho=7850.0, cp=480.0)

result = bf.AdjointCGM(data, material=mat).solve()
```
