## FZG LOAD STAGES ############################################################
def gtorque(load, arm):
    if arm == '0.35':
        if load == 'k01':
            return 3.3
        elif load == 'k02':
            return 13.7
        elif load == 'k03':
            return 28.785
        elif load == 'k04':
            return 46.635     
        elif load == 'k05':
            return 69.98
        elif load == 'k06':
            return 98.82    
        elif load == 'k07':
            return 132.455
        elif load == 'k08':
            return 171.585
        elif load == 'k09':
            return 215.513
        elif load == 'k10':
            return 265
        elif load == 'k11':
            return 319.18
        elif load == 'k12':
            return 378.26
        elif load == 'k13':
            return 438.8533
        elif load == 'k14':
            return 499.94
    elif arm == '0.50':
        if load == 'k01':
            return 3.3
        elif load == 'k02':
            return 13.7
        elif load == 'k03':
            return 35.25
        elif load == 'k04':
            return 60.75     
        elif load == 'k05':
            return 94.1
        elif load == 'k06':
            return 135.3    
        elif load == 'k07':
            return 183.35
        elif load == 'k08':
            return 239.25
        elif load == 'k09':
            return 302
        elif load == 'k10':
            return 372.7
        elif load == 'k11':
            return 450.1
        elif load == 'k12':
            return 534.5
        elif load == 'k13':
            return 626.9333
        elif load == 'k14':
            return 714.2