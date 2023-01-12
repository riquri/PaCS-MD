#!/bin/sh

name=initial

acpype -p ./leap.prmtop -x ./leap.inpcrd


#375  Na+  NA+ 5862   6.554   1.802   6.141
#sed -i -e "s/  Na+  NA+/  Na+   IP/g" ./leap_GMX.gro
#sed -i -e "s/  Cl-  Cl-/  Cl-   IM/g" ./leap_GMX.gro


#     1       IP      1          NA+         NA+       1      1     22.9898
#     1       IP      1          NA+          IP       1      1     22.9898
#sed -i -e "s/IP      1          NA+         NA+/IP      1          NA+          IP/g" ./leap_GMX.top
#sed -i -e "s/IM      1          CL-         CL-/IM      1          CL-          IM/g" ./leap_GMX.top


# Na+      Na+         0.00000  0.00000   A     2.43928e-01   3.65846e-01 ; 1.37  0.0874
#sed -i -e "s/Na+      Na+/IP       IP/g" ./leap_GMX.top
#sed -i -e "s/Cl-      Cl-/IM       IM/g" ./leap_GMX.top


cp ./leap.amb2gmx/leap_GMX.top ./${name}.top
cp ./leap.amb2gmx/leap_GMX.gro ./${name}.gro



