# a = length in relative units of the ccd
# b = length in relative units of the "dead gap"

# "GAP FROM HANDBOOK FIG 4.1"
a = 2048*(15e-6)  # 15 micron square pixels
b = 1.28e-3 # Figure 4.1 gap

A_tot = 4*a**2 + 4*a*b + b**2
A_dead = 4*a*b + b**2

print('Model with 1.28 mm gap handbook Figure 4.1')
print(f'A_dead/A_tot = {A_dead/A_tot*100:.2f}%')
print(f'1-A_dead/A_tot = {(1-A_dead/A_tot)*100:.2f}%')

# "GAP FROM HANDBOOK TEXT"
a = 2048*(15e-6)  # 15 micron square pixels
b = 1.5e-3 # Figure 4.1 gap

A_tot = 4*a**2 + 4*a*b + b**2
A_dead = 4*a*b + b**2

print('Model with 1.5 mm gap handbook Figure 4.1')
print(f'A_dead/A_tot = {A_dead/A_tot*100:.2f}%')
print(f'1-A_dead/A_tot = {(1-A_dead/A_tot)*100:.2f}%')


# "SIMPLE MODEL"
a = 2048
b = (2/0.015)

A_tot = 4*a**2 + 4*a*b + b**2
A_dead = 4*a*b + b**2

print('Simple Model')
print(f'A_dead/A_tot = {A_dead/A_tot*100:.2f}%')
print(f'1-A_dead/A_tot = {(1-A_dead/A_tot)*100:.2f}%')
