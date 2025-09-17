import argparse
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def add_filter_column_change(dataframe):
    
    """
    Add a filter change column to the df based on the following conditions:
    The new column 'FILTER_change' will indicate:
    - 'No change PASS' if both FILTER_truth and FILTER_query contain 'PASS'
    - 'No change EXCLUDE' if both FILTER_truth and FILTER_query are 'EXCLUDE'
    - 'Change from PASS to EXCLUDE' if FILTER_truth contains 'PASS'
       and FILTER_query contains 'EXCLUDE'
    - 'Change from EXCLUDE to PASS' if FILTER_truth contains 'EXCLUDE'
       and FILTER_query contains 'PASS'
    - 'Variant removed PASS' if FILTER_truth contains 'PASS'
       and FILTER_query is '.'
    - 'Variant removed EXCLUDE' if FILTER_truth contains 'EXCLUDE'
       and FILTER_query is '.'
    - 'Variant added PASS' if FILTER_truth is '.'
       and FILTER_query contains 'PASS'
    - 'Variant added EXCLUDE' if FILTER_truth is '.'
       and FILTER_query contains 'EXCLUDE'
    - 'No change, rescued variant' if either FILTER truth or query contains
      'rescued' and are identical
    - 'Changed rescued variant' if either FILTER truth or query contains
      'rescued' and are not identical
    
    Args:
        dataframe (pd.DataFrame): The input dataframe containing 'FILTER_truth'
                                  and 'FILTER_query' columns.
    
    Returns:
        pd.DataFrame: The modified dataframe with the new 'FILTER_change'
                      column.
    """
    # Set up all conditions
    conditions = [
        (dataframe['FILTER_truth'].str.contains('PASS')) &
        (dataframe['FILTER_query'].str.contains('PASS')),
        
        (dataframe['FILTER_truth'].str.contains('EXCLUDE')) &
        (dataframe['FILTER_query'].str.contains('EXCLUDE')),
        
        (dataframe['FILTER_truth'].str.contains('PASS')) &
        (dataframe['FILTER_query'].str.contains('EXCLUDE')),
        
        (dataframe['FILTER_truth'].str.contains('EXCLUDE')) &
        (dataframe['FILTER_query'].str.contains('PASS')),
        
        (dataframe['FILTER_truth'].str.contains('PASS')) & 
        (dataframe['FILTER_query'] == '.'),
        
        (dataframe['FILTER_truth'].str.contains('EXCLUDE')) &
        (dataframe['FILTER_query'] == '.'),
        
        (dataframe['FILTER_truth'] == '.') &
        (dataframe['FILTER_query'].str.contains('PASS')),
        
        (dataframe['FILTER_truth'] == '.') &
        (dataframe['FILTER_query'].str.contains('EXCLUDE')),
        
        (dataframe['FILTER_truth'].str.contains('rescued') |
         dataframe['FILTER_query'].str.contains('rescued')) & 
        (dataframe['FILTER_truth'] == dataframe['FILTER_query']),
        
        (dataframe['FILTER_truth'].str.contains('rescued') |
         dataframe['FILTER_query'].str.contains('rescued')) & 
        (dataframe['FILTER_truth'] != dataframe['FILTER_query'])
    ]
    
    # Assign strings for each condition specified
    choices = [
        'No change PASS',
        'No change EXCLUDE',
        'Change from PASS to EXCLUDE',
        'Change from EXCLUDE to PASS',
        'Variant removed PASS',
        'Variant removed EXCLUDE',
        'Variant added PASS',
        'Variant added EXCLUDE',
        'No change, rescued variant',
        'Change rescued variant'
    ]
    
    # Create new column with the new conditions
    dataframe['FILTER_change'] = np.select(conditions, choices,
                                           default='Unknown')
    
    return dataframe


def create_interactive_correlation_plot(csv_file):
    """
    Create an interactive correlation plot using VAF_truth vs VAF_query 
    as x and y axes and colored by the FILTER_change
    
    Parameters
    ----------
    csv_file : str
        Filepath corresponding to the sample.merged.csv
    
    Returns
    -------
    fig : plotly.graph_objs._figure.Figure
        Plotly figure of the created correlation plot
    """
    try:
        # Read the CSV file
        df = pd.read_csv(csv_file)
        df = add_filter_column_change(df)

        print("The unique values for FILTER_truth",
              df['FILTER_truth'].unique())
        print("The unique values for FILTER_query",
              df['FILTER_query'].unique())
        
        # Check if required columns exist
        required_cols = ['VAF_truth', 'VAF_query', 'QUAL.truth', 'CHROM']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"Error: Missing columns: {missing_cols}")
            print(f"Available columns: {list(df.columns)}")
            return
        
        
        # Calculate correlation
        correlation = df['VAF_truth'].corr(df['VAF_query'])
        print("Checking customdata order:")
        print("hover_data keys:", list(df.columns))
        
        # Define custom colour mapping
        color_discrete_map = {
            'No change PASS': "#00CC96",
            'No change EXCLUDE': "#636EFA",
            'Change from PASS to EXCLUDE': "#FFA15A",
            'Change from EXCLUDE to PASS': "#FECB52",
            'Variant removed PASS': "#F14526",
            'Variant removed EXCLUDE': "#19D3F3",
            'Variant added PASS': '#FF6692',
            'Variant added EXCLUDE': "#AB63FA",
            'No change, rescued variant': "#00CC1B",
            'Change rescued variant': "#A1A31F",
            'Unknown': '#696969'
        }
        
        # Create interactive scatter plot
        sample = df['sample'][0]
        fig = px.scatter(df,
                         x='VAF_truth',
                         y='VAF_query',
                         color='FILTER_change',
                         color_discrete_map=color_discrete_map,
                         hover_data={
                             'CHROM': True,
                             'POS': True,
                             'REF': True,
                             'ALT': True,
                             'VAF_truth': True,
                             'VAF_query': True,
                             'DP_truth': True,
                             'DP_query': True,
                             'FILTER_truth': True,
                             'FILTER_query': True,
                             'FILTER_change': True
                         },
                         color_continuous_scale='viridis',
                         title=f'VAF Correlation Plot for {sample}<br>Pearson '
                               f'Correlation: {correlation:.3f}',
                         labels={
                            'FILTER_change': 'Filter Change'
                         }
                         )
        print("Checking the first row customdata:",
              fig.data[0].customdata[0] 
              if fig.data[0].customdata is not None else "None")
        
        # Update text box when hovering over datapoints
        fig.update_traces(
            hovertemplate='<b>Variant: %{customdata[0]}:%{customdata[1]} ' +
                          '%{customdata[2]}>%{customdata[3]}</b><br>' +
                          '<b>VAF_truth:</b> %{x}<br>' +
                          '<b>VAF_query:</b> %{y}<br>' +
                          '<b>DP_truth:</b> %{customdata[4]}<br>' +
                          '<b>DP_query:</b> %{customdata[5]}<br>' +
                          '<b>FILTER_truth:</b> %{customdata[6]}<br>' +
                          '<b>FILTER_query:</b> %{customdata[7]}<br>' +
                          '<b>FILTER_change:</b> %{customdata[8]}<br>' +
                          '<extra></extra>'
        )
        
        # Add diagonal reference line
        min_val = min(df['VAF_truth'].min(), df['VAF_query'].min())
        max_val = max(df['VAF_truth'].max(), df['VAF_query'].max())
        
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
        
        # Add grid on the plot output
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray',
                         title='VAF Flowcell truth')
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray',
                         title='VAF Flowcell query')
        
        return fig
    
    # Raise an error if the processing has not taken place
    except Exception as e:
        print(f"Error processing {csv_file}: {str(e)}")
        return None


def save_plot_pdf_html(input_file, output):
    """
    Create correlation plots for all CSV files and save to PDF using Plotly
    Parameters
    ----------
    input_file : str
        Filepath corresponding to the sample.merged.csv
    output : str
        Filepath corresponding to the output .html and .pdf files
    """
    fig = create_interactive_correlation_plot(input_file)
    
    # Save as HTML
    fig.write_html(f"{output}.html")
    print(f"Interactive plot saved to {output}.html")
    
    # Save as pdf
    fig.write_image(f"{output}.pdf")
    print(f"Interactive plot saved to {output}.pdf")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate .html and .pdf from the sample.merged.csv files"
    )
    parser.add_argument("--feature_file", help="Input feature .merged.csv file")
    parser.add_argument("--output", help="Output filepath, no need to specify .pdf or .html")
    args = parser.parse_args()
    
    return args


def main():
    """
    Main function
    """
    args = parse_args()
    save_plot_pdf_html(args.feature_file, args.output)


if __name__ == "__main__":
    main()
