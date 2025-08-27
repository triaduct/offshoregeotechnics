from shiny import App, render, ui
import auxiliary.averaging_methods
import pandas as pd, numpy as np
import matplotlib.pyplot as plt

def give_qcavg_profile(averaging_method,diameter,z,qc):
    cpt = pd.DataFrame({"z":z,"qc":qc})
    cpt["qcavg"] = np.nan
    
    cpt.z.apply(lambda z: auxiliary.averaging_methods.lcpc(cpt,diameter,z))
    
    

app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_select("averaging_method", "Averaging method", choices = ["Koppejan 4D/8D", "LCPC 1.5D", "Boulanger & DeJong","De Boorder"]),
            ui.input_slider("D", "Diameter", min=5, max=30, value=10, step=1),
        ),
        ui.output_plot("plot"),
    ),
)


def server(input, output, session):
    @output
    @render.plot(alt="Sine wave")
    
    def plot():
        x,y = give_qcavg_profile(input.averaging_method(), input.D())

        fig, ax = plt.subplots()
        ax.plot(x,y)
        ax.set_xlim(-1,1.5)
        ax.set_ylim(0,1.5)
        return fig

app = App(app_ui, server)

