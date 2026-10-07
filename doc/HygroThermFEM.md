```mermaid
graph TD
    %% Package definitions
    HygroThermFEM[HygroThermFEM<br/><font color='red'>v_1.0.75</font>]
    KeffCavity[KeffCavity<br/><font color='red'>Version_1.1.8</font>]
    Eigen[Eigen<br/><font color='red'>5.0.1</font>]
    WindowsCalcEngine[Windows-CalcEngine<br/><font color='red'>Version_1.0.77</font>]

    %% Package relations
    HygroThermFEM --> KeffCavity
    HygroThermFEM --> Eigen
    KeffCavity --> WindowsCalcEngine
```
