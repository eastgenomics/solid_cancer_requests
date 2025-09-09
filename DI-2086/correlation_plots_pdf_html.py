import argparse
import glob
import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

def add_filter_column_change(dataframe):
    """
    Add a filter change column to the dataframe based on the conditions specified.
    The new column 'FILTER_change' will indicate:
    - 'No change PASS' if both FILTER_s1 and FILTER_s2 contain 'PASS'
    - 'No change EXCLUDE' if both FILTER_s1 and FILTER_s2 contain 'EXCLUDE'
    - 'Change from PASS to EXCLUDE' if FILTER_s1 contains 'PASS' and FILTER_s2 contains 'EXCLUDE'
    - 'Change from EXCLUDE to PASS' if FILTER_s1 contains 'EXCLUDE' and FILTER_s2 contains 'PASS'
    - 'Variant removed PASS' if FILTER_s1 contains 'PASS' and FILTER_s2 is '.'
    - 'Variant removed EXCLUDE' if FILTER_s1 contains 'EXCLUDE' and FILTER_s2 is '.'
    - 'Variant added PASS' if FILTER_s1 is '.' and FILTER_s2 contains 'PASS'
    - 'Variant added EXCLUDE' if FILTER_s1 is '.' and FILTER_s2 contains 'EXCLUDE'

    Args:
        dataframe (pd.DataFrame): The input dataframe containing 'FILTER_s1' and 'FILTER_s2' columns.

    Returns:
        pd.DataFrame: The modified dataframe with the new 'FILTER_change' column.
    """
    conditions = [
        (dataframe['FILTER_s1'].str.contains('PASS')) & (dataframe['FILTER_s2'].str.contains('PASS')),
        (dataframe['FILTER_s1'].str.contains('EXCLUDE')) & (dataframe['FILTER_s2'].str.contains('EXCLUDE')),
        (dataframe['FILTER_s1'].str.contains('PASS')) & (dataframe['FILTER_s2'].str.contains('EXCLUDE')),
        (dataframe['FILTER_s1'].str.contains('EXCLUDE')) & (dataframe['FILTER_s2'].str.contains('PASS')),
        (dataframe['FILTER_s1'].str.contains('PASS')) & (dataframe['FILTER_s2'] == '.'),
        (dataframe['FILTER_s1'].str.contains('EXCLUDE')) & (dataframe['FILTER_s2'] == '.'),
        (dataframe['FILTER_s1'] == '.') & (dataframe['FILTER_s2'].str.contains('PASS')),
        (dataframe['FILTER_s1'] == '.') & (dataframe['FILTER_s2'].str.contains('EXCLUDE')),
    ]
    choices = [
        'No change PASS',
        'No change EXCLUDE',
        'Change from PASS to EXCLUDE',
        'Change from EXCLUDE to PASS',
        'Variant removed PASS',
        'Variant removed EXCLUDE',
        'Variant added PASS',
        'Variant added EXCLUDE'
    ]
    
    dataframe['FILTER_change'] = np.select(conditions, choices, default='Unknown')
    
    return dataframe

def create_interactive_correlation_plot(csv_file):
    """
    Create an interactive correlation plot with VAF_s1 vs VAF_s2 colored by tag
    with CHROM information displayed on hover
    #TO BE CHANGED
    
    Todo's
    - Print out the unique values for the filter column to verify possible combination
    - Write a function that adds a new column based on the conditions specified
    - Make it plot and colour appropriately
    - Present work
    """
    try:
        # Read the CSV file
        df = pd.read_csv(csv_file)
        
        df = add_filter_column_change(df)
        
        print("Columns of the dataframe: ", df.columns)
        print("The unique values for FILTER_s1", df['FILTER_s1'].unique())
        print("The unique values for FILTER_s2", df['FILTER_s2'].unique())
        
        # Check if required columns exist
        required_cols = ['VAF_s1', 'VAF_s2', 'QUAL.truth', 'CHROM']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"Error: Missing columns: {missing_cols}")
            print(f"Available columns: {list(df.columns)}")
            return
        
        
        # Calculate correlation
        correlation = df['VAF_s1'].corr(df['VAF_s2'])
        print("Checking customdata order:")
        print("hover_data keys:", list(df.columns))
        
        # Define custom colour mapping
        color_discrete_map = {
            'No change PASS': "#00CC96",
            'No change EXCLUDE': "#636EFA",
            'Change from PASS to EXCLUDE': "#FFA15A",
            'Change from EXCLUDE to PASS': "#FECB52",
            'Variant removed PASS': "#EF553B",
            'Variant removed EXCLUDE': "#19D3F3",
            'Variant added PASS': '#FF6692',
            'Variant added EXCLUDE': "#AB63FA",
            'Unknown': '#696969'
        }
        
        # Create interactive scatter plot
        sample = df['sample'][0]
        fig = px.scatter(df,
                         x='VAF_s1',
                         y='VAF_s2',
                         color='FILTER_change',
                         color_discrete_map=color_discrete_map,
                         hover_data={
                             'CHROM': True,
                             'POS': True,
                             'REF': True,
                             'ALT': True,
                             'VAF_s1': True,
                             'VAF_s2': True,
                             'DP_s1': True,
                             'DP_s2': True,
                             'FILTER_s1': True,
                             'FILTER_s2': True,
                             'FILTER_change': True
                         },
                         color_continuous_scale='viridis',
                         title=f'VAF Correlation Plot for {sample}<br>Pearson Correlation: {correlation:.3f}',
                         labels={
                            'FILTER_change': 'Filter Change'
                         }
                         )
        print("First row customdata:", fig.data[0].customdata[0] if fig.data[0].customdata is not None else "None")
        fig.update_traces(
            hovertemplate='<b>Variant: %{customdata[0]}:%{customdata[1]} %{customdata[2]}>%{customdata[3]}</b><br>' +
                          '<b>VAF_s1:</b> %{x}<br>' +
                          '<b>VAF_s2:</b> %{y}<br>' +
                          '<b>DP_s1:</b> %{customdata[4]}<br>' +
                          '<b>DP_s2:</b> %{customdata[5]}<br>' +
                          '<b>FILTER_s1:</b> %{customdata[6]}<br>' +
                          '<b>FILTER_s2:</b> %{customdata[7]}<br>' +
                          '<b>FILTER_change:</b> %{customdata[8]}<br>' +
                          '<extra></extra>'
        )
        
        # Add diagonal reference line
        min_val = min(df['VAF_s1'].min(), df['VAF_s2'].min())
        max_val = max(df['VAF_s1'].max(), df['VAF_s2'].max())
        
        fig.for_each_trace(
            lambda trace: trace.update(name=trace.name.replace("TP", "Shared")
                                       .replace("FP", "S2 only")
                                       .replace("FN", "S1 only"))
            if trace.name in ["TP", "FP", "FN"] else trace
)
        
        fig.add_trace(go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode='lines',
            line=dict(dash='dash', color='red', width=1),
            name='Pearson Correlation line',
            hovertemplate='Reference line y=x<extra></extra>'
        ))
        
        # Update layout
        fig.update_layout(
            width=1200,
            height=900,
            showlegend=True,
            font=dict(size=12)
        )
        
        # Add grid
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', title='VAF Flowcell S1')
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', title='VAF Flowcell S2')
        
        return fig
        
    except Exception as e:
        print(f"Error processing {csv_file}: {str(e)}")
        return None

def create_plot_pdf_html(input_file, output):
    """
    Create correlation plots for all CSV files and save to PDF using Plotly
    """
    fig = create_interactive_correlation_plot(input_file)
    
    # Save as HTML
    fig.write_html(f"{output}.html")
    print(f"Interactive plot saved to {output}.html")
    
    # Save as pdf
    fig.write_image(f"{output}.pdf")
    print(f"Interactive plot saved to {output}.pdf")


def main(features_file, output):
    """
    Main function
    """
    create_plot_pdf_html(features_file, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate .html and .pdf from the sample.merged.csv files"
    )
    parser.add_argument("--feature_file", help="Input feature .merged.csv file")
    parser.add_argument("--output", help="Output filepath, no need to specify .pdf or .html")
    args = parser.parse_args()
    main(args.feature_file, args.output)
    