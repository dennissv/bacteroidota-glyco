#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 15:10:06 2024

@author: dennis

Identify which monosaccharides are present
    Diagnostic ions
    Peptide delta-delta masses

155.0 - Addition to Prevotella species

"""

from itertools import combinations_with_replacement
import urllib.parse

import pandas as pd
import numpy as np
import streamlit as st
from glycowork.motif.draw import GlycoDraw

PROTON_MASS = 1.007276
PPM_ERROR = 50
# MONOSACCHARIDES = ('HexNAc', 'Hex', 'dHex')
ANNOTATIONS_IONS = {'HexNAc-C2H6O3': 126.055, # https://www.nature.com/articles/s43586-022-00128-4/tables/4
                    'HexNAc-CH6O3': 138.055,
                    'HexNAc-C2H4O2': 144.0655,
                    'HexNAc-2H2O': 168.0655,
                    'HexNAc-H2O': 186.0761,
                    # 'HexNAc(+)': 204.0866, # Modified according to diagnostic ion mining paper
                    'Hex-2H2O': 127.0390,
                    'Hex-H2O': 145.0495,
                    # 'Hex(+)': 163.0601,
                    'HexP': 243.0270,
                    'NeuAc(+)': 292.1027,
                    'NeuGc(+)': 308.0976,
                    }

ANNOTATIONS_PEPTIDES = {'dHex': 146.057909,
                        # 'dHexMe': 146.057909 + 14.015650,
                        'Hex': 162.052824,
                        # 'HexP': 162.052824 + 79.966331,
                        'HexA': 176.03209,
                        'HexNAc': 203.079373,
                        # 'Unknown2': 192.023,
                        # 'Unknown': 182.10019999999997,
                        # 'Unknown': 202.1074-PROTON_MASS,
                        # 'Unknown': 155.0,
                        # 'Unknown': 258.1084 - PROTON_MASS,
                        # 'Unknown': 171.40074933333335,
                        # 'NeuAc': 291.095417,
                        # 'NeuGc': 308.0976 - PROTON_MASS,
                        # 'NeuGc': 42.01, # Acetylation
                        # 'Fuc': 14.015650, # Methylation
               }

def within_ppm(v1, v2, ppm=50):
    if abs(v1 - v2) <= (ppm * v1 / 1e6):
        return True
    return False

def ppm_error(v1, v2):
    return round(abs(v1 - v2) / v1 * 1e6, 2)

def merge_values_within_ppm(values, ppm=50):
    groups = dict()
    for v1 in values:
        found = False
        for v2 in groups:
            if within_ppm(v1, v2, ppm=ppm):
                groups[v2].append(v1)
                found = True
                break
        if not found:
            groups[v1] = [v1]
    
    out = sorted([(len(x), np.average(x)) for x in groups.values()], key=lambda x: x[0], reverse=True)
    
    return out

# def read_diagnostic_ions(file_path):
#     lines = []
#     with open(file_path) as f:
#         for line in f:
#             lines.append(line.strip().split('\t'))
#     data = pd.DataFrame(lines[1:], columns=lines[0], index=None)
#     float_columns = list(data.columns[9:21]) + ['Intensity'] + list(data.columns[29:32])
#     data[float_columns]= data[float_columns].apply(pd.to_numeric, errors='coerce')
#     filtered_df = data[data['Peptide'].str.contains('DT|DS')]
#     return filtered_df

def find_monosaccharides(df):
    monosaccharides = dict()
    monosaccharides_delta_masses = dict()
    for index, row in df.iterrows():
        peak = row['mass']
        ion_type = row['ion_type']
        # if ion_type == 'peptide':
        #     for name, mass in ANNOTATIONS_PEPTIDES.items():
        #         if within_ppm(peak, mass):
        #             df.at[index, 'annotation'] = name
        #             df.at[index, 'ppm_error'] = ppm_error(peak, mass)
                    
        if ion_type == 'diagnostic':
            for name, mass in ANNOTATIONS_IONS.items():
                if within_ppm(peak, mass):
                    # df.at[index, 'annotation'] = name
                    # df.at[index, 'ppm_error'] = ppm_error(peak, mass)
                    monosaccharides[name] = mass
    
    peptide_masses = df[df['ion_type'] == 'peptide']
    peptide_masses = sorted(list(peptide_masses['mass']) + [row['peak_apex']], reverse=True)
    delta_masses = []
    for i, m1 in enumerate(peptide_masses):
        for m2 in peptide_masses[i+1:]:
            delta_mass = m1 - m2
            delta_masses.append(delta_mass)
    for m1 in delta_masses:
        for name, m2 in ANNOTATIONS_PEPTIDES.items():
            if within_ppm(m1, m2, PPM_ERROR):
                if name not in monosaccharides_delta_masses:
                    monosaccharides_delta_masses[name] = 0
                monosaccharides_delta_masses[name] += 1
    # st.write(delta_masses)
    for value in merge_values_within_ppm(delta_masses, ppm=300):
        st.write(value)
    # st.write(monosaccharides_delta_masses)
    
    return monosaccharides

def delta_delta_masses(df):
    pass

df = None
peptide_masses = []
st.title('GlycAnalyser')

# img = GlycoDraw('Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Neu5Gc(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc', show_linkage=False)
# svg_string = img.as_svg()
# st.markdown(f"<div style='text-align: center;'>{svg_string}</div>", unsafe_allow_html=True)

uploaded_file = st.sidebar.file_uploader("Choose a file")

if uploaded_file is not None:
    df = pd.read_csv(uploaded_file, delimiter='\t')
    df = df[df['mass'] >= 120]
    df = df.drop(columns=['mod_annotation', 'remainder_propensity', 'avg_charge', 'intensity_fold_change'])
    df['annotation'] = ''
    df['ppm_error'] = ''
    
    unique_peaks = pd.unique(df['peak_apex'])
    for peak in unique_peaks:
        new_row = df.iloc[-1].copy()
        new_row['peak_apex'] = peak
        new_row['mass'] = peak
        new_row['ion_type'] = 'peptide'
        new_row['percent_mod'] = 100
        new_row['percent_unmod'] = 0
        new_row['avg_intensity_mod'] = 100
        new_row['delta_mod_mass'] = 0
        
        new_row_df = pd.DataFrame([new_row])
        df = pd.concat([df, new_row_df], ignore_index=True)
    
    diagnostic_ion_combinations = dict()
    for r in range(1, 8):  # From 1 element up to 4 elements
        for combo in combinations_with_replacement(ANNOTATIONS_PEPTIDES.keys(), r):
            mass = sum(ANNOTATIONS_PEPTIDES[monosaccharide] for monosaccharide in combo) + PROTON_MASS
            diagnostic_ion_combinations['-'.join(combo) + '(+)'] = mass
    
    peptide_ion_combinations = dict()
    for r in range(1, 12):  # From 1 element up to 4 elements
        for combo in combinations_with_replacement(ANNOTATIONS_PEPTIDES.keys(), r):
            mass = sum(ANNOTATIONS_PEPTIDES[monosaccharide] for monosaccharide in combo)
            peptide_ion_combinations['-'.join(combo)] = mass
    
    for index, row in df.iterrows():
        peak = row['mass']
        ion_type = row['ion_type']
        best = np.inf
        if ion_type == 'peptide':
            for name, mass in peptide_ion_combinations.items():
                current_ppm_error = ppm_error(peak, mass)
                if (current_ppm_error < PPM_ERROR) and (current_ppm_error < best):
                    df.at[index, 'annotation'] = name
                    df.at[index, 'ppm_error'] = current_ppm_error
                    best = current_ppm_error
                    
        elif ion_type == 'diagnostic':
            for name, mass in {**diagnostic_ion_combinations, **ANNOTATIONS_IONS}.items():
                current_ppm_error = ppm_error(peak, mass)
                if (current_ppm_error < 20) and (current_ppm_error < best):
                    df.at[index, 'annotation'] = name
                    df.at[index, 'ppm_error'] = ppm_error(peak, mass)
    df.sort_values(by=['peak_apex', 'mass'], ascending=[False, True], inplace=True)
    cols = df.columns.tolist()
    new_order = cols[:2] + ['annotation', 'ppm_error'] + cols[2:-2]
    df = df[new_order]

# Sidebar
if type(df) != type(None):
    option = st.sidebar.selectbox('Select a peak:', (x for x in pd.unique(df['peak_apex']) if x > 300))
    selected_rows = st.sidebar.checkbox('Show only selected peak')

# Filtering data based on sidebar selection
    if selected_rows:
        filtered_df = df[df['peak_apex'] == option]
        peptide_masses = df[df['peak_apex'] == option]
    else:
        filtered_df = df
    filtered_df = filtered_df.reset_index(drop=True)

    monosaccharides = find_monosaccharides(filtered_df)
    # for name, mass in monosaccharides.items():
    #     st.write(f'{name}, {mass}')
    
    svg_strings = []
    for index, row in filtered_df.iterrows():
        if row['ion_type'] == 'diagnostic':
            if '(+)' in row['annotation']:
                glyco_string = '(1-6)'.join(row['annotation'][:-3].split('-'))
                img = GlycoDraw(glyco_string, show_linkage=False, compact=True)
                svg = img.as_svg()
                svg_string = f"data:image/svg+xml;utf8,{urllib.parse.quote(svg)}"
                svg_strings.append(svg_string)
                # st.write(glyco_string)
            else:
                svg = f'<svg width="350" height="100" xmlns="http://www.w3.org/2000/svg"><text x="10" y="50" font-family="Arial" font-size="40" fill="black">{row["annotation"]}</text></svg>'
                svg_string = f"data:image/svg+xml;utf8,{urllib.parse.quote(svg)}"
                svg_strings.append(svg_string)
        else:
            if row['annotation']:
                glyco_string = '(1-6)'.join(row['annotation'].split('-'))
                img = GlycoDraw(glyco_string, show_linkage=False, compact=True)
                svg = img.as_svg()
                svg_string = f"data:image/svg+xml;utf8,{urllib.parse.quote(svg)}"
                svg_strings.append(svg_string)
            else:
                svg_strings.append('')
    # encoded_svgs = [f"data:image/svg+xml;utf8,{urllib.parse.quote(svg)}" for svg in svg_strings]
    filtered_df.insert(3, 'SNFG', svg_strings)
    filtered_df = filtered_df.drop(columns=('annotation'))
    # st.write(svg_strings)

    st.data_editor(
    data=filtered_df,
    column_config={
        'SNFG': st.column_config.ImageColumn(label='SNFG', width="medium")
    },
    hide_index=True
    )
    
    # peptide_ions = filtered_df[filtered_df['ion_type'] == 'peptide']
    # peptide_ions = sorted(list(peptide_ions['mass']), reverse=True)
    # delta_masses
    # for i, m1 in enumerate(peptide_ions):
    #     for m2 in peptide_ions[i+1:]:
            
    # for index, row in peptide_ions.iterrows():
    #     st.write(index)
    
    # Displaying the DataFrame
    # st.write(filtered_df)
    # st.data_editor(
    #     data=filtered_data_figures,
    #     column_config={
    #         "SVG Images": st.column_config.ImageColumn(label="SVG Preview", width="medium")
    #     },
    #     hide_index=True
    # )
