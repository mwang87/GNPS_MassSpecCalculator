# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table
import plotly.express as px
from dash.dependencies import Input, Output
import os
from zipfile import ZipFile
import urllib.parse
from flask import Flask
import json
import pandas as pd
from molmass import Formula
import requests
import urllib.parse


server = Flask(__name__)
app = dash.Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server


NAVBAR = dbc.Navbar(
    children=[
        dbc.NavbarBrand(
            html.Img(src="https://gnps-cytoscape.ucsd.edu/static/img/GNPS_logo.png", width="120px"),
            href="https://gnps.ucsd.edu"
        ),
        dbc.Nav(
            [
                dbc.NavItem(dbc.NavLink("GNPS Mass Spec Calculator", href="#")),
            ],
        navbar=True)
    ],
    color="light",
    dark=False,
    sticky="top",
)

DASHBOARD = [
    dbc.CardHeader(html.H5("GNPS Mass Spec Calculator")),
    dbc.CardBody(
        [
            html.Div(id='version', children="Version - Release_1"),
            html.Br(),
            dbc.Label("Molecular Formula", html_for="formula_entry"),
            dbc.Input(className="mb-3", id='formula_entry', placeholder="Enter Formula"),
            html.Br(),
            dbc.Label("Smiles Structure", html_for="smiles_entry"),
            dbc.Textarea(className="mb-3", id='smiles_entry', placeholder="Enter SMILES structure"),
            html.Hr(),
            html.Div(children="Adduct Table"),
            html.Br(),
            dcc.Loading(
                className="mb-3",
                id="massspecinfo",
                children=[html.Div([html.Div(id="loading-output-3")])],
                type="default",
            )
        ]
    )
]

BODY = dbc.Container(
    [
        dbc.Row([dbc.Col(dbc.Card(DASHBOARD)),], style={"marginTop": 30}),
    ],
    className="mt-12",
)

app.layout = html.Div(children=[NAVBAR, BODY])

# This function will rerun at any 
@app.callback(
    [Output('massspecinfo', 'children')],
    [Input('formula_entry', 'value'), Input('smiles_entry', 'value')],
)
def generate_url(formula_entry, smiles_entry):
    exact_mass = 0

    if formula_entry is not None and len(formula_entry):
        f = Formula(formula_entry)
        exact_mass = f.isotope.mass
    else:
        # Getting exact mass
        url = "https://gnps-structure.ucsd.edu/structuremass?smiles={}".format(urllib.parse.quote(smiles_entry))
        r = requests.get(url)
        exact_mass = float(r.text)

    adducts_to_report = ["M", "M+H", "M-H"]
    output_list = []

    for adduct in adducts_to_report:
        adduct_mass = get_adduct_mass(exact_mass, adduct)
        output_dict = {}
        output_dict["adduct"] = adduct
        output_dict["mz"] = adduct_mass
        output_list.append(output_dict)

    table_fig = dash_table.DataTable(
        columns=[
            {"name": i, "id": i, "deletable": True, "selectable": True} for i in ["adduct", "mz"]
        ],
        data=output_list,
        editable=True,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        column_selectable="single",
        selected_columns=[],
        selected_rows=[],
        page_action="native",
        page_current= 0,
        page_size= 10,
    )

    return [table_fig]



def get_adduct_mass(exact_mass, adduct):
    M = exact_mass

    if adduct == 'M':
        return M, 1

    if adduct == ('M+3H'):
        return M/3 + 1.007276, 3

    if adduct == ('M+2H+Na'):
        return M/3 + 8.334590, 3

    if adduct == ('M+H+2Na'):
        return M/3 + 15.7661904, 3

    if adduct == ('M+3Na'):
        return M/3 + 22.989218, 3

    if adduct == ('M+2H'):
        return M/2 + 1.007276, 2

    if adduct == ('M+H+NH4'):
        return M/2 + 9.520550, 2

    if adduct == ('M+H+Na'):
        return M/2 + 11.998247, 2

    if adduct == ('M+H+K'):
        return M/2 + 19.985217, 2

    if adduct == ('M+ACN+2H'):
        return M/2 + 21.520550, 2

    if adduct == ('M+2Na'):
        return M/2 + 22.989218, 2

    if adduct == ('M+2ACN+2H'):
        return M/2 + 42.033823, 2

    if adduct == ('M+3ACN+2H'):
        return M/2 + 62.547097, 2

    if adduct == ('M+H'):
        return M + 1.007276, 1

    if adduct == ('M+H-H2O'):
        return M + 19.01839 + 1.007276 + 1.007276, 1

    if adduct == ('M+NH4'):
        return M + 18.033823, 1

    if adduct == ('M+Na'):
        return M + 22.989218, 1

    if adduct == ('M+CH3OH+H'):
        return M + 33.033489, 1

    if adduct == ('M+K'):
        return M + 38.963158, 1

    if adduct == ('M+ACN+H'):
        return M + 42.033823, 1

    if adduct == ('M+2Na-H'):
        return M + 44.971160, 1

    if adduct == ('M+IsoProp+H'):
        return M + 61.06534, 1

    if adduct == ('M+ACN+Na'):
        return M + 64.015765, 1

    if adduct == ('M+2K-H'):
        return M + 76.919040, 1

    if adduct == ('M+DMSO+H'):
        return M + 79.02122, 1

    if adduct == ('M+2ACN+H'):
        return M + 83.060370, 1

    if adduct == ('M+IsoProp+Na+H'):
        return M + 84.05511, 1

    if adduct == ('2M+H'):
        return 2*M + 1.007276, 1

    if adduct == ('2M+NH4'):
        return 2*M + 18.033823, 1

    if adduct == ('2M+Na'):
        return 2*M + 22.989218, 1

    if adduct == ('2M+K'):
        return 2*M + 38.963158, 1

    if adduct == ('2M+ACN+H'):
        return 2*M + 42.033823, 1

    if adduct == ('2M+ACN+Na'):
        return 2*M + 64.015765, 1

    if adduct == ('M-H2O+H'):
        return M - 17.00384, 1

    if adduct == ('M-3H'):
        return M/3 - 1.007276, -3

    if adduct == ('M-2H'):
        return M/2 - 1.007276, -2

    if adduct == ('M-H2O-H'):
        return M - 19.01839, -1

    if adduct == ('M-H'):
        return M - 1.007276, -1

    if adduct == ('M+Na-2H'):
        return M + 20.974666, -2

    if adduct == ('M+Cl'):
        return M + 34.969402, 1

    if adduct == ('M+K-2H'):
        return M + 36.948606, -2

    if adduct == ('M+FA-H'):
        return M + 44.998201, -1

    if adduct == ('M+Hac-H'):
        return M + 59.013851, -1

    if adduct == ('M+Br'):
        return M + 78.918885, 1

    if adduct == ('M+TFA-H'):
        return M + 112.985586, -1

    if adduct == ('2M-H'):
        return 2*M - 1.007276, -1

    if adduct == ('2M+FA-H'):
        return 2*M + 44.998201, -1

    if adduct == ('2M+Hac-H'):
        return 2*M + 59.013851, -1

    if adduct == ('3M-H'):
        return 3*M - 1.007276, -1

    if adduct == ('M-2H2O+H'):
        return M + 1.007276 - 2*18.01057, 1

    if adduct == ('2M-2H+Na'):
        return M*2 - 1.007276 *2 + 22.989218, -1
    if adduct == ('2M-2H+K'):
        return M*2 - 1.007276 *2 + 38.963158, -1

    return exact_mass, 0

if __name__ == "__main__":
    app.run_server(debug=True, port=5000, host="0.0.0.0")
