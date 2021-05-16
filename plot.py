import matplotlib.pyplot as p
import numpy as n
import matplotlib as m
m.rcParams.update({'font.size':15})

y = lambda x : 0.025*x + 0.25


x = n.linspace(-10, 0, 50)
p.plot(x, y(x), color="red", linewidth=2.0, label="other GCs calibration line")
p.plot(-2.9, 0.2857, color="blue", markersize=20.0, marker="*", label="Horologium I")
p.ylim(0, 0.5)
p.xlim(-9,-2)
p.xlabel("Absolute GC magnitude")
p.ylabel("Fraction of binaries")
p.legend()
p.tight_layout()
p.show()