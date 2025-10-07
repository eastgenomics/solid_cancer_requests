import argparse

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from sklearn.metrics import r2_score


def add_filter_column_change(merged_var_df):
    """
    Add a filter change column to the dataframe based on the following conditions:
    The new column 'FILTER_change' will indicate:
    - 'No change PASS' if both FILTER_truth and FILTER_query contain 'PASS'
        Example pair values: 
            - 'PASS' -> 'PASS'
            - 'rescued' -> 'PASS'
            - 'PASS' -> 'rescued'
            - 'LowSupport;LowDP;rescued' -> 'LowSupport;LowDP;rescued'
            - 'rescued' -> 'rescued'
    - 'No change EXCLUDE' if both FILTER_truth and FILTER_query have 'EXCLUDE'
        Example pair values:
            - 'EXCLUDE' -> 'EXCLUDE'
            - 'LowSupport;LowDP;rescued;EXCLUDE' -> 'LowSupport;LowDP;rescued;EXCLUDE'
            - 'LowSupport;EXCLUDE' -> 'EXCLUDE'
            - 'EXCLUDE' -> 'LowSupport;EXCLUDE'
    - 'Change from PASS to EXCLUDE' if FILTER_truth contains 'PASS'
       and FILTER_query contains 'EXCLUDE'
        Example pair values:
            - 'PASS' -> 'EXCLUDE'
            - 'PASS' -> 'LowSupport;EXCLUDE'
            - 'PASS' -> 'LowSupport;LowDP;rescued;EXCLUDE'
            - 'rescued' -> 'EXCLUDE'
            - 'LowSupport;rescued' -> 'LowSupport;LowDP;rescued;EXCLUDE'
    - 'Change from EXCLUDE to PASS' if FILTER_truth contains 'EXCLUDE'
       and FILTER_query contains 'PASS'
        Example pair values:
            - 'EXCLUDE' -> 'PASS'
            - 'LowSupport;EXCLUDE' -> 'PASS'
            - 'LowSupport;LowDP;rescued;EXCLUDE' -> 'PASS'
            - 'EXCLUDE' -> 'rescued'
            - 'LowSupport;EXCLUDE' -> 'LowSupport;rescued'
    - 'Variant removed PASS' if FILTER_truth contains 'PASS'
       and FILTER_query is '.'
        Example pair values:
            - 'PASS' -> '.'
            - 'rescued' -> '.'
    - 'Variant removed EXCLUDE' if FILTER_truth contains 'EXCLUDE'
       and FILTER_query is '.'
        Example pair values:
            - 'EXCLUDE' -> '.'
            - 'LowSupport;EXCLUDE' -> '.'
            - 'LowSupport;LowDP;rescued;EXCLUDE' -> '.'
    - 'Variant added PASS' if FILTER_truth is '.'
       and FILTER_query contains 'PASS'
        Example pair values:
            - '.' -> 'PASS'
    - 'Variant added EXCLUDE' if FILTER_truth is '.'
       and FILTER_query contains 'EXCLUDE'
        Example pair values:
            - '.' -> 'EXCLUDE'
            - '.' -> 'LowSupport;EXCLUDE'
            - '.' -> 'LowSupport;LowDP;rescued;EXCLUDE'

    NOTES:
        - If more than one condition is met, the first condition in the list
          will be applied.
        - If none of the conditions are met, 'Unknown' will be assigned.
        - No cases where FILTER could be NaN or 'rescued;PASS' have been
          observed in the datasets, so these have not been explicitly handled.
    Args:
        merged_var_df (pd.DataFrame): The input dataframe containing 
                                      'FILTER_truth' and 'FILTER_query'
                                      columns.

    Returns:
        merged_var_df (pd.DataFrame): The modified dataframe with the new 'FILTER_change'
                                           column.
    """
    # Set up all conditions
    conditions = [
        # No change PASS
        (
            merged_var_df["FILTER_truth"].str.contains("PASS|rescued") 
            & ~merged_var_df["FILTER_truth"].str.contains("EXCLUDE")
        )
        & (
            merged_var_df["FILTER_query"].str.contains("PASS|rescued")
            & ~merged_var_df["FILTER_query"].str.contains("EXCLUDE")
        ),
        # No change EXCLUDE
        (merged_var_df["FILTER_truth"].str.contains("EXCLUDE"))
        & (merged_var_df["FILTER_query"].str.contains("EXCLUDE")),
        # Change from PASS to EXCLUDE
        (
            merged_var_df["FILTER_truth"].str.contains("PASS|rescued") 
            & ~merged_var_df["FILTER_truth"].str.contains("EXCLUDE")
        )
        & (merged_var_df["FILTER_query"].str.contains("EXCLUDE")),
        # Change from EXCLUDE to PASS
        (merged_var_df["FILTER_truth"].str.contains("EXCLUDE"))
        & (
            merged_var_df["FILTER_query"].str.contains("PASS|rescued")
            & ~merged_var_df["FILTER_query"].str.contains("EXCLUDE")
        ),
        # Variant removed PASS
        (
            merged_var_df["FILTER_truth"].str.contains("PASS|rescued") 
            & ~merged_var_df["FILTER_truth"].str.contains("EXCLUDE")
        )
        & (merged_var_df["FILTER_query"] == "."),
        # Variant removed EXCLUDE
        (merged_var_df["FILTER_truth"].str.contains("EXCLUDE"))
        & (merged_var_df["FILTER_query"] == "."),
        # Variant added PASS
        (merged_var_df["FILTER_truth"] == ".")
        & (
            merged_var_df["FILTER_query"].str.contains("PASS|rescued")
            & ~merged_var_df["FILTER_query"].str.contains("EXCLUDE")
        ),
        # Variant added EXCLUDE
        (merged_var_df["FILTER_truth"] == ".")
        & (merged_var_df["FILTER_query"].str.contains("EXCLUDE")),
    ]

    # Assign strings for each condition specified
    choices = [
        "No change PASS",
        "No change EXCLUDE",
        "Change from PASS to EXCLUDE",
        "Change from EXCLUDE to PASS",
        "Variant removed PASS",
        "Variant removed EXCLUDE",
        "Variant added PASS",
        "Variant added EXCLUDE"
    ]

    # Create new column with the new conditions
    merged_var_df["FILTER_change"] = np.select(conditions, choices,
                                               default="Unknown")

    return merged_var_df


def create_interactive_correlation_plot(csv_file, output_path):
    """
    Create an interactive correlation plot using VAF_truth vs VAF_query
    as x and y axes and colored by the FILTER_change

    Parameters
    ----------
    csv_file : str
        Filepath corresponding to the sample.merged.csv
    output_path : str
        Filepath corresponding to the output_path for the .html and .pdf files

    Returns
    -------
    fig : plotly.graph_objs._figure.Figure
        Plotly figure of the created correlation plot
    """
    try:
        # Read the CSV file
        records_df = pd.read_csv(csv_file)

        # Check if required columns exist
        required_cols = [
            "CHROM", "POS", "REF", "REF.truth", "ALT", "ALT.truth",
            "DP_truth", "DP_query", "FILTER_query", "FILTER_truth",
            "VAF_truth", "VAF_query", "sample"
        ]
        missing_cols = [
            col for col in required_cols if col not in records_df.columns
            ]

        if missing_cols:
            print(f"Available columns: {list(records_df.columns)}")
            raise ValueError(f"Missing columns: {missing_cols}")

        df = add_filter_column_change(records_df)

        print("The unique values for FILTER_truth", df["FILTER_truth"].unique())
        print("The unique values for FILTER_query", df["FILTER_query"].unique())

        # Calculate correlation of all samples
        r2_total = r2_score(df["VAF_truth"], df["VAF_query"])
        print(f"R² correlation: {r2_total:.4f}")
        
        # Calculate the correlation of only PASS related samples
        subset_df = df[df["FILTER_change"].isin(
            ["No change PASS",
             "Change from PASS to EXCLUDE",
             "Change from EXCLUDE to PASS",
             "Variant added PASS",
             "Variant removed PASS"])]
        r2_subset = r2_score(subset_df["VAF_truth"], subset_df["VAF_query"])
        print(f"R² correlation (PASS related): {r2_subset:.4f}")
        print("Checking customdata order:")
        print("hover_data keys:", list(df.columns))

        # Define custom colour mapping
        color_discrete_map = {
            "No change PASS": "#00CC96",
            "No change EXCLUDE": "#636EFA",
            "Change from PASS to EXCLUDE": "#FFA15A",
            "Change from EXCLUDE to PASS": "#FECB52",
            "Variant removed PASS": "#F14526",
            "Variant removed EXCLUDE": "#19D3F3",
            "Variant added PASS": "#FF6692",
            "Variant added EXCLUDE": "#AB63FA",
            "Unknown": "#696969",
        }

        # Create interactive scatter plot
        sample = df["sample"].iloc[0]
        fig = px.scatter(
            df,
            x="VAF_truth",
            y="VAF_query",
            color="FILTER_change",
            color_discrete_map=color_discrete_map,
            custom_data=['CHROM','POS','REF','ALT',
                         'DP_truth','DP_query',
                         'FILTER_truth','FILTER_query',
                         'FILTER_change'],
            category_orders={
                "FILTER_change": list(color_discrete_map.keys())
            },
            title=f"VAF Correlation Plot for {sample}<br>"
            f"R² Correlation of all samples: {r2_total:.4f} <br>"
            f"R² Correlation of PASS related samples: {r2_subset:.4f}",
            labels={"FILTER_change": "Filter Change"},
        )
        print(
            "Checking the first row customdata:",
            fig.data[0].customdata[0] if fig.data[0].customdata is not None else "None",
        )

        # Update text box when hovering over datapoints
        fig.update_traces(
            hovertemplate="<b>Variant: %{customdata[0]}:%{customdata[1]} "
            + "%{customdata[2]}>%{customdata[3]}</b><br>"
            + "<b>VAF_truth:</b> %{x}<br>"
            + "<b>VAF_query:</b> %{y}<br>"
            + "<b>DP_truth:</b> %{customdata[4]}<br>"
            + "<b>DP_query:</b> %{customdata[5]}<br>"
            + "<b>FILTER_truth:</b> %{customdata[6]}<br>"
            + "<b>FILTER_query:</b> %{customdata[7]}<br>"
            + "<b>FILTER_change:</b> %{customdata[8]}<br>"
            + "<extra></extra>"
        )

        # Add diagonal reference line
        min_val = min(df["VAF_truth"].min(), df["VAF_query"].min())
        max_val = max(df["VAF_truth"].max(), df["VAF_query"].max())

        fig.add_trace(
            go.Scatter(
                x=[min_val, max_val],
                y=[min_val, max_val],
                mode="lines",
                line=dict(dash="dash", color="red", width=1),
                name="y=x reference line",
                hovertemplate="Reference line y=x<extra></extra>",
            )
        )

        # Update layout
        fig.update_layout(width=1200, height=900, showlegend=True, font=dict(size=12))

        # Add grid on the plot output
        fig.update_xaxes(
            showgrid=True, gridwidth=1, gridcolor="lightgray", title="VAF truth"
        )
        fig.update_yaxes(
            showgrid=True, gridwidth=1, gridcolor="lightgray", title="VAF query"
        )
        # Output plot df for fine-resolution inspection
        df.to_csv(f"{output_path}_with_FILTER_change.merged.csv")
        return fig

    # Raise an error if the processing has not taken place
    except Exception as e:
        print(f"Error processing {csv_file}: {e!s}")
        raise


def save_plot_pdf_html(input_file, output_path):
    """
    Create correlation plots for all CSV files and save to PDF using Plotly
    Parameters
    ----------
    input_file : str
        Filepath corresponding to the sample.merged.csv
    output_path : str
        Filepath corresponding to the output .html and .pdf files
    """
    fig = create_interactive_correlation_plot(input_file, output_path)
    if fig is None:
        raise RuntimeError(f"Failed to create figure from {input_file}")

    # Save as HTML
    fig.write_html(f"{output_path}.html")
    print(f"Interactive plot saved to {output_path}.html")

    # Save as pdf
    try:
        fig.write_image(f"{output_path}.pdf")
        print(f"Static pdf saved to {output_path}.pdf")
    except Exception as e:
        raise RuntimeError(
            "Failed to export PDF. Ensure 'kaleido' is installed. "
            f"Error: {e}"
        ) from e



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate .html and .pdf from the sample.merged.csv files"
    )
    parser.add_argument("--feature_file", required=True,
                        help="Input feature .merged.csv file")
    parser.add_argument("--output_path", required=True,
                        help="Output filepath, no need to specify "
                             ".pdf or .html")
    args = parser.parse_args()

    return args


def main():
    """
    Main function
    """
    args = parse_args()
    save_plot_pdf_html(args.feature_file, args.output_path)


if __name__ == "__main__":
    main()
