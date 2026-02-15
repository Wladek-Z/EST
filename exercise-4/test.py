import numpy as np
import pandas as pd

if __name__ == "__main__":
    v = pd.read_csv("Vinet_params.txt")['V0'].values
    v_b = v[0] / 0.529177249**3
    a = (4*v_b)**(1/3)
    ar = a / np.sqrt(2)
    print(ar)