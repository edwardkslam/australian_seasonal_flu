import matplotlib as mpl
import matplotlib.pyplot as plt

plt.style.use("seaborn-white")
mpl.rcParams['axes.labelsize']= "x-large"
mpl.rcParams['xtick.labelsize']= "x-large"
mpl.rcParams['ytick.labelsize']= "x-large"
mpl.rcParams['axes.titlesize']= "x-large"
mpl.rcParams['axes.formatter.use_mathtext'] = True
mpl.rcParams['axes.formatter.limits'] = ((-3, 3))
mpl.rcParams['text.usetex'] = False
mpl.rcParams['font.sans-serif'] = ["Arial",
                                   "Helvetica",
                                   "Liberation Sans"]
mpl.rcParams['font.family'] = "sans-serif"

lw = 5

letter_loc = (-0.1, 1.15)
letter_size = "x-large"
