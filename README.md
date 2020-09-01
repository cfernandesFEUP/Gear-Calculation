The software is capable to calculate gear geometry according to MAAG book, predict Hertz contact pressure, gear and rolling bearings power losses - only NJ 406 MA and QJ 308 N2MA are implemented because the software is intended for usage with FZG test rig results!

Other Features not included on the repository: gear CAD geometry generation and automatic structured FEM mesh, CalculiX themo-mechanical integration, heat transfer coefficients calculation - next release.

This version includes:

- the main file is 'GearCP.py':

    'GearCP.py' allows to select:
    
        - Gear geometry (C14, C40, H501, H701, H951) which is stored into the file "gears.py" - additional geometries can be added
        
        - Gear material:
        
            mat = [STEEL STEEL]; meaning pinion and gear material respectively;
            
        - Operating Conditions:
        
            nmotor - list of motor speeds (meaning output speed)
            load - list of FZG load stages, available with load arm of 0.35 m and 0.5 m (check LoadStage.py); a required torque can also be defined.
            
        - Oil Selection:
        
           'dry' - no lubricant, a Cofficient of Friction should be given - useful for plastic gearing
           'PAOR', 'MINR', etc - "oils.py"  according to papers [1, 2, 3]. It calculates the Coefficient of Friction using Schlenk approach and corresponding XL -  additional lubricants can be added.
           
           
 How to Cite?

 This program is on the basis of the following papers:
 
 [1] Fernandes, C. M. C. G., Martins, R. C., & Seabra, J. H. O. (2014). Torque loss of type C40 FZG gears lubricated with wind turbine gear oils. Tribology International, 70(0), 83–93. https://doi.org/10.1016/j.triboint.2013.10.003
 
 [2] Fernandes, C. M. C. G., Marques, P. M. T., Martins, R. C., & Seabra, J. H. O. (2015). Gearbox power loss. Part I: Losses in rolling bearings. Tribology International, 88(0), 298–308. https://doi.org/10.1016/j.triboint.2014.11.017
 
 [3] Fernandes, C. M. C. G., Marques, P. M. T. T., Martins, R. C., & Seabra, J. H. O. (2015). Gearbox power loss. Part II: Friction losses in gears. Tribology International, 88, 309–316. https://doi.org/10.1016/j.triboint.2014.12.004
 
 [4] Fernandes, C. M. C. G., Marques, P. M. T., Martins, R. C., & Seabra, J. H. O. (2015). Gearbox power loss. Part III: Application to a parallel axis and a planetary gearbox. Tribology International, 88, 317–326. https://doi.org/10.1016/j.triboint.2015.03.029
 
 [5] Fernandes, C., Martins, R., Seabra, J. H. O., & Blazquez, L. (2016). FZG Gearboxes Lubricated with Different Formulations of Polyalphaolefin Wind Turbine Gear Oils. August, 56–60.
 
 [6] Fernandes, C. M. C. G. M., Hammami, M., Martins, R. C., & Seabra, J. H. O. H. (2016). Power loss prediction: Application to a 2.5 MW wind turbine gearbox. Proceedings of the Institution of Mechanical Engineers, Part J: Journal of Engineering Tribology, 230(8), 983–995. https://doi.org/10.1177/1350650115622362
 
 Dry conditions/polymer gears:
 
 [7] Fernandes, C. M. C. G., Rocha, D. M. P., Martins, R. C., Magalhães, L., & Seabra, J. H. O. (2018). Finite element method model to predict bulk and flash temperatures on polymer gears. Tribology International, 120, 255–268. https://doi.org/10.1016/j.triboint.2017.12.027
 
 [8] Fernandes, C. M. C. G., Rocha, D. M. P., Martins, R. C., Magalhães, L., & Seabra, J. H. O. (2019). Hybrid Polymer Gear Concepts to Improve Thermal Behavior. Journal of Tribology, 141(3), 032201. https://doi.org/10.1115/1.4041461
