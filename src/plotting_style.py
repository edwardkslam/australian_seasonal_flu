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
mpl.rcParams['font.sans-serif'] = ["Helvetica",
                                   "Arial",
                                   "Liberation Sans"]
mpl.rcParams['font.family'] = "sans-serif"

lw = 5
