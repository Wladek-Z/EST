B_0 = float(input("Enter B_0: "))
mols = int(input("Enter number of molecules: "))
B_0_convert = round(B_0 * 1e-9 * 1e30 * (13.6056931229947 * 1.60218e-19) / mols, 6)
print(f"{B_0_convert} GPa")
