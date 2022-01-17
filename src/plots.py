from pathlib import Path
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import os
from ipywidgets import (Tab, SelectMultiple, Accordion, ToggleButton,
                        VBox, HBox, HTML, Button, Image, Output)
from ipywidgets import HBox, VBox, Layout, HTML, Dropdown
from src.cif_plots import JsMolFigure
from IPython.display import display
#from IPython.display import Image


class SolubilityPlot:

    def __init__(self, csv_path: str) -> None:
        self.df: pd.DataFrame = pd.read_csv(csv_path)

    def change_amine(self, change):
        if change['name'] == 'label':
            self.update_plot(change['new'])

    def get_data(self, amine):
        df = self.df[self.df['Amine'] == amine]
        data = []
        for i, row in df.sort_values(axis=0, by='Solvent').iterrows():
            data.append(go.Scatter(x=[row['Max PbI [M]']],
                                   y=[row['Max Amine [M]']],
                                   name=row['Solvent'],
                                   mode='markers',
                                   marker=dict(size=12),
                                   )
                        )
        return data

    def update_plot(self, amine):
        sol_data = self.get_data(amine)
        with self.sol_fig.batch_update():
            for i, trace in enumerate(sol_data):
                self.sol_fig.data[i].x = trace.x
                self.sol_fig.data[i].y = trace.y

    def plot(self):
        amines = self.df['Amine'].unique()
        traces = self.get_data(amines[0])
        options = [(amine, i) for i, amine in enumerate(amines)]
        self.dropdown_amines = Dropdown(
            options=options, value=0, description='Amine: ')
        self.dropdown_amines.observe(self.change_amine)
        self.sol_fig = go.FigureWidget(data=traces)
        self.sol_fig.update_layout(
            title="Solubility Plot",
            xaxis_title="Max PbI [M]",
            yaxis_title="Max Amine [M]",
            legend_title="Solvents",
            font=dict(size=18),
            width=1000,
            height=600,
        )
        return VBox([self.dropdown_amines, self.sol_fig])


class ReactionOutcomeWidget:

    def __init__(self, reaction_outcome_path: str):
        self.df = pd.read_csv(reaction_outcome_path,
                              error_bad_lines=False, skipinitialspace=True)
        self.df.fillna('', inplace=True)
        self.reaction_summary = pd.read_csv('./data/reaction_summary.csv')

    def get_data(self, rxn_id, ):
        #diffusion_heights_folder = Path('./data/Diffusion_rate_and_crystal_height_all_CSV')
        #filename = list(diffusion_heights_folder.glob(f'get_liquid_height_and_crystal_start_{rxn_id}*.*'))
        csv_folder = Path('./data/csv_for_notebook')
        solvent_ht_filepath = csv_folder / f'{rxn_id}_solvent_height.csv'
        conc_filepath = csv_folder / f'{rxn_id}_concentrations.csv'

        sol_data = []
        conc_data = []
        max_sol = 0
        max_conc = 0
        if solvent_ht_filepath.exists():
            #filename = filename[0]
            sol_df = pd.read_csv(solvent_ht_filepath)
            conc_df = pd.read_csv(conc_filepath)
            antisol = sol_df['AS Height Experimental (cm)']
            solvent = sol_df['Solvent Height Experimental (cm)']
            sol_data = [go.Scatter(x=sol_df['Time for Height Build Up (hours)'],
                                   y=solvent, name='ML', mode='markers+lines',),
                        go.Scatter(x=sol_df['Time for Height Build Up (hours)'],
                                   y=antisol, name='AS', mode='markers+lines'),
                        ]
            max_sol = np.max(np.concatenate([solvent, antisol]))
            conc_data = [go.Scatter(x=conc_df['Time for Concentration (hours)'], y=conc_df['Solvent [M]'], name='S', mode='markers+lines'),
                         go.Scatter(x=conc_df['Time for Concentration (hours)'],
                                    y=conc_df['DCM [M]'], name='AS', mode='markers+lines'),
                         go.Scatter(x=conc_df['Time for Concentration (hours)'],
                                    y=conc_df['FAH [M]'], name='FAH', mode='markers+lines'),
                         go.Scatter(x=conc_df['Time for Concentration (hours)'],
                                    y=conc_df['Inorganic [M]'], name='Inorganic', mode='markers+lines'),
                         go.Scatter(x=conc_df['Time for Concentration (hours)'], y=conc_df['Organic [M]'], name='Organic', mode='markers+lines'), ]
            max_conc = np.max(np.concatenate(
                [conc_df['Solvent [M]'], conc_df['DCM [M]']]))
            protonation_time = sol_df.iloc[0]['Crystallization Time (hours)']

        else:
            print(f'Missing file for {rxn_id}')
        return sol_data, conc_data, max_sol, max_conc, protonation_time

    def update_fig(self, rxn_id, ):
        sol_data, conc_data, max_sol, max_conc, protonation_time = self.get_data(
            rxn_id, )
        with self.sol_fig.batch_update():
            for i, trace in enumerate(sol_data):
                self.sol_fig.data[i].x = trace.x
                self.sol_fig.data[i].y = trace.y
        self.sol_fig.update_shapes(
            {'x0': protonation_time, 'y0': 0, 'x1': protonation_time, 'y1': max_sol, })
        with self.conc_fig.batch_update():
            for i, trace in enumerate(conc_data):
                self.conc_fig.data[i].x = trace.x
                self.conc_fig.data[i].y = trace.y

        self.conc_fig.update_shapes(
            {'x0': protonation_time, 'y0': 0, 'x1': protonation_time, 'y1': max_conc, })
        base_image_path = Path(
            f'./data/visible_nucleation_images/{rxn_id}.png')
        if base_image_path.exists():
            with base_image_path.open('rb') as f:
                image = f.read()
                #self.nuc_image = Image(value=image, format='png', width=600, height=400)
                self.nuc_image.value = image

    def render_figs(self, rxn_id):
        sol_data, conc_data, max_sol, max_conc, protonation_time = self.get_data(
            rxn_id)
        layout = go.Layout(title='Test')
        self.sol_fig = go.FigureWidget(data=sol_data, layout=layout)
        self.sol_fig.update_layout(
            title="Solvent Height",
            xaxis_title="Time (hours)",
            yaxis_title="Meniscus height (cm)",
            # legend_title="Solvents",
            font=dict(size=18),
            width=900,
            height=700,
        )
        self.sol_fig.add_shape(type="line",
                               x0=protonation_time, y0=0, x1=protonation_time, y1=max_sol,
                               line=dict(color="RoyalBlue", width=3)
                               )
        self.conc_fig = go.FigureWidget(data=conc_data)
        self.conc_fig.update_layout(
            title="Concentration Plot",
            xaxis_title="Time (hrs)",
            yaxis_title='Concentration [M]',
            font=dict(size=18),
            width=900,
            height=700,
        )
        self.conc_fig.add_shape(type="line",
                                x0=protonation_time, y0=0, x1=protonation_time, y1=max_conc,
                                line=dict(color="RoyalBlue", width=3)
                                )
        with open("./data/visible_nucleation_images/MA_338_2.png", "rb") as f:
            image = f.read()
            self.nuc_image = Image(
                value=image, format='png', width=600, height=400)

        #self.out = Output(layout={'border': '1px solid black'})
        #self.jsmol = JsMolFigure(['./data/cifs/ajn19-019.cif'], ['./data/cifs/ajn19-019.cif'], {}, widget_side=600)

        # with self.out:
        #    display(self.jsmol.plot)

    def render_plot(self):
        self.button_grid = self.generate_button_grid()
        self.render_figs('MA_354_1')
        # return VBox([self.button_grid, HBox([self.sol_fig, self.conc_fig,]), HBox([self.nuc_image])])
        self.items = [self.sol_fig, self.conc_fig, self.nuc_image]
        self.tab = Tab()
        self.tab.children = self.items
        self.tab._titles = {0: 'Solvent Height Plot',
                            1: 'Concentration Plot', 2: 'First visible nucleation'}

        with open('./data/rapid2_TABLE.PNG', "rb") as f:
            image = f.read()
            self.table = Image(value=image, format='png', width=350)

        return VBox([HBox([self.button_grid, self.table]), self.tab])

    def on_button_clicked(self, b):
        #print(f"Button clicked. {b.amine} {b.solvent}")
        rxn_id = self.reaction_summary[(self.reaction_summary['Amine'] == b.amine)
                                       & (self.reaction_summary['Solvent'] == b.solvent)]['Rxn_ID']
        if not rxn_id.empty:
            rxn_id = rxn_id.iloc[0]
            self.sol_fig.update_layout(
                title=f"Solubility Plot for {b.amine} {b.solvent}",)
            self.update_fig(rxn_id,)

    def generate_button_grid(self):
        hbox_list = []
        for col in self.df.columns:
            vbox = [HTML(f'<b>{col}</b>')]
            for i, row in enumerate(self.df[col]):
                description = ''
                button_style = ''
                if row:
                    if len(row.split(',')) == 1:
                        element = HTML(f'<b>{row}</b>')
                    elif len(row.split(',')) == 2:
                        description, button_style = row.split(',')
                        element = Button(
                            description=description,
                            disabled=False,
                            button_style=button_style.strip(),
                            layout=Layout(width='50px'),
                        )
                        element.on_click(self.on_button_clicked)
                        element.amine = col
                        element.solvent = self.df['Solvent'].iloc[i]
                vbox.append(element)
            hbox_list.append(VBox(vbox, layout=Layout(border='1px solid')))
        return HBox(hbox_list)


"""
traces = []
        self.layout = go.Layout(
            scene=dict(
                xaxis=dict(
                    title='Solvent',
                    tickmode='linear',
                    dtick=0.5,
                    #range=[0, self.max_inorg],
                ),
                yaxis=dict(
                    title='Max Amine [M]',
                    tickmode='linear',
                    dtick=0.5,
                    #range=[0, self.max_org],
                ),
                zaxis=dict(
                    title='Max PbI [M]',
                    tickmode='linear',
                    dtick=1.0,
                    #range=[0, self.max_acid],
                ),
            ),
            legend=go.layout.Legend(
                x=0,
                y=1,
                traceorder="normal",
                font=dict(
                    family="sans-serif",
                    size=12,
                    color="black"
                ),
                bgcolor="LightSteelBlue",
                bordercolor="Black",
                borderwidth=2
            ),
            width=975,
            height=700,
            margin=go.layout.Margin(
                l=10,
                r=10,
                b=10,
                t=10,
                pad=1
            ),
        )
        for amine in amines:
            filtered_df = self.df[self.df['Amine']==amine]
            solvent = filtered_df['Solvent']
            max_amine = filtered_df['Max Amine [M]']
            max_pbi = filtered_df['Max PbI [M]']
            trace = go.Scatter3d(
                    x=solvent,
                    y=max_amine,
                    z=max_pbi,
                    mode='markers',
                    name=f'{amine}',
                    text=[f"<b>Amine</b>: {row['Amine']}"
                          f"<br><b>Solvent</b>: {row['Solvent']}," 
                          f"<br><b>Max Amine [M]</b>: {row['Max Amine [M]']},"
                          f"<br><b>Max PbI [M]</b>: {row['Max PbI [M]']}"
                        for i, row in filtered_df.iterrows()
                        ],
                    hoverinfo='text',
                    marker=dict(
                        size=4,
                        #color=self.trace_colors[i],
                        line=dict(
                            width=0.2
                        ),
                        opacity=1.0
                    )
                )
            traces.append(trace)
"""
