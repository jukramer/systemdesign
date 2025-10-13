import math

g = 9.80665
R = 287
T0 = 288.15
p0 = 101325
ro = 1.225
h0 = 0

layers = [
    (11000, -0.0065),
    (20000, 0),
    (32000, 0.001),
    (47000, 0.0028),
    (51000, 0),
    (71000, -0.0028),
    (86000, -0.002),]

print('\n     **** ISA calculator ****\n')
print("1. Altitude in meters\n2. Altitude in feet\n3. Altitude in FL")

while True:
    pick = int(input("Enter your choice: "))
    
    if pick == 1:
        alt = float(input("\nEnter altitude [m]: "))
        while True:
                if alt > 86000:
                    print("Sorry , I can only do altitudes up to 86000 m.")
                    alt = float(input("Altitude [m]: "))
                else:
                    break
        break
    elif pick == 2:
        alt = float(input("\nEnter altitude [ft]: "))
        while True:
                if alt > 282152.23:
                    print("Sorry , I can only do altitudes up to 282152.23 ft.")
                    alt = float(input("Altitude [ft]: "))
                else:
                    alt = alt/3.281
                    break
        break
    elif pick == 3:
        alt = float(input("\nEnter altitude [FL]: "))
        while True:
                if alt > 2821.5223:
                    print("Sorry , I can only do altitudes up to 2821.5223 FL.")
                    alt = float(input("Altitude [FL]: "))
                else:
                    alt = alt/3.281*100
                    break
        break
    else:
        print("Please choose option 1, 2 or 3")

print("\nWould you like to use a different sea level temperature? (if not, press enter)")
T_user = input("Temperature: ")
print(T_user)

if T_user != "":
    T0 = float(T_user) + 273.15

for hLimit, a in layers:
    if alt <= h0:
        break

    h1 = min(alt, hLimit)
    delta_h = h1-h0

    if a == 0:
        T1 = T0
        p1 = p0 * math.exp((-g*delta_h/(R*T0)))
    else:
        T1 = T0 + a * delta_h
        p1 = p0 * (T1/T0)**(-g/(a*R))

    T0 = T1
    p0 = p1
    h0 = h1
    
density = p0/(R*T0)

print("Temperature:", round(T0, 2), "K (", round(T0-273.15, 2), "'C)")
print("Pressure:", round(p0), "Pa (", round(p0/101325*100), "% SL)" )
print("Density:", round(density, 4), "kg/m3 (", round(density/ro*100), "% SL)")

dummy = input("\nPress enter to end the ISA calculator.")
