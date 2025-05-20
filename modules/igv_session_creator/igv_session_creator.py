import marimo

__generated_with = "0.12.0"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
async def _():
    try:
        import micropip
        await micropip.install("jinja2")
    except Exception:
        ...
    return (micropip,)


@app.cell(hide_code=True)
def _(mo):
    from io import StringIO
    from jinja2 import Environment, FileSystemLoader
    import pandas as pd
    import time
    import asyncio

    # Initialize jinja2 environment
    template_folder = mo.notebook_location() / "public"
    env = Environment(loader=FileSystemLoader(searchpath=f"{template_folder}"))
    session_template = env.get_template("igv_session_template.xml.j2")
    tableau10_rgb = [
        "31,119,180",   # Blue
        "255,127,14",   # Orange
        "44,160,44",    # Green
        "214,39,40",    # Red
        "148,103,189",  # Purple
        "140,86,75",    # Brown
        "227,119,194",  # Pink
        "127,127,127",  # Gray
        "188,189,34",   # Olive
        "23,190,207"    # Cyan
    ]
    tableau10_rgb = {
        "Blue": "31,119,180",
        "Orange": "255,127,14",
        "Green": "44,160,44",
        "Red": "214,39,40",
        "Purple": "148,103,189",
        "Brown": "140,86,75",
        "Pink": "227,119,194",
        "Gray": "127,127,127",
        "Olive": "188,189,34",
        "Cyan": "23,190,207"
    }

    tableau10_rgb_form = {_name:  mo.md(f"<div style='display:inline-block;width:40px;height:20px;background-color:rgb({_rgb});border:1px solid #ccc'></div>") for _i, (_name, _rgb) in enumerate(tableau10_rgb.items()) }
    return (
        Environment,
        FileSystemLoader,
        StringIO,
        asyncio,
        env,
        pd,
        session_template,
        tableau10_rgb,
        tableau10_rgb_form,
        template_folder,
        time,
    )


@app.cell(hide_code=True)
def _(mo):
    # Load dataframe in from file in TSV format
    fh = mo.ui.file(filetypes=[".tsv", ".tsv.gz"], multiple=True, kind="area")
    mo.vstack([mo.ui.file(kind="button"), fh])
    return (fh,)


@app.cell(hide_code=True)
def _(StringIO, fh, mo, pd):
    mo.stop(not fh.value)

    decoded_str = fh.contents(0).decode("utf-8")
    metadata = pd.read_csv(StringIO(decoded_str), sep="\t")
    return decoded_str, metadata


@app.cell(hide_code=True)
def _(mo, pd):
    tracks = pd.DataFrame(columns=["name", "path", "type", "color"])
    get_tracks, set_tracks = mo.state(tracks)
    return get_tracks, set_tracks, tracks


@app.cell(hide_code=True)
def _(mo):
    genome = mo.ui.text(placeholder="mm39", label="Genome")
    genome
    return (genome,)


@app.cell(hide_code=True)
def _(genome, metadata, mo):
    metadata_filtered = mo.ui.dataframe(metadata)
    mo.vstack([mo.md(f"::lucide:dna::Selected Genome: **{genome.value}**"), metadata_filtered])
    return (metadata_filtered,)


@app.cell(hide_code=True)
def _(mo, tableau10_rgb_form):
    color_selection = mo.ui.table(
        data=tableau10_rgb_form,
        initial_selection=[0],
        selection="single",
        show_download=False,
    )
    return (color_selection,)


@app.cell(hide_code=True)
def _(color_selection, metadata, mo, tableau10_rgb_form):
    columns = metadata.columns.tolist()

    selected_color = color_selection.value[0]["key"]

    name_col = mo.ui.dropdown(columns, label="::lucide:pencil:: **Name column**")
    url_col = mo.ui.dropdown(columns, label="::lucide:link:: **URL column**")
    color = mo.hstack([mo.md(f"::lucide:palette:: **Track Color**:"), tableau10_rgb_form[selected_color]], justify="start")
    type_col = mo.ui.dropdown(["bw"], label="::lucide:chart-no-axes-column:: **Track type column**")
    button = mo.ui.run_button(label="Add Tracks")

    mo.hstack([mo.vstack([name_col, url_col, color, type_col, button], justify="space-between"),color_selection])
    return button, color, columns, name_col, selected_color, type_col, url_col


@app.cell(hide_code=True)
def _(pd):
    def add_new_tracks(df, new_rows):
        return pd.concat([df, new_rows], ignore_index=True)

    def remove_rows_by_index(df, indices_to_remove):
        return df.drop(index=indices_to_remove).reset_index(drop=True)


    #def add_tracks(df, rows_to_add):
    #    return (
    #        pd.concat([df, rows_to_add], ignore_index=True)
    #        .drop_duplicates(subset=["name"])
    #        .assign(source="metadata")
    #    )
    return add_new_tracks, remove_rows_by_index


@app.cell(hide_code=True)
def _(
    add_new_tracks,
    button,
    metadata_filtered,
    mo,
    name_col,
    pd,
    selected_color,
    set_tracks,
    tableau10_rgb,
    type_col,
    url_col,
):
    mo.stop(not button.value)

    new_rows = pd.DataFrame.from_dict([{
                "name": row[name_col.value],
                "path": row[url_col.value],
                "type": type_col.value,
                "color": tableau10_rgb[selected_color],
            }for _, row in metadata_filtered.value.iterrows() ])

    set_tracks(lambda df: df.pipe(add_new_tracks, new_rows=new_rows))
    return (new_rows,)


@app.cell(hide_code=True)
def _(current_view, delete_button, mo, remove_rows_by_index, set_tracks):
    mo.stop(not delete_button.value or current_view.value.empty)

    set_tracks(lambda df: df.pipe(remove_rows_by_index, indices_to_remove=current_view.value.index))
    return


@app.cell
def _(mo):
    mo.md(
        """
        Add "clear all"
        use dataframe widget
        add autoscale
        track height
        name of the track - suffix/prefix
        multi-table support
        terra integration
        """
    )
    return


@app.cell
def _(get_tracks, mo):
    current_view = mo.ui.table(get_tracks())
    delete_button = mo.ui.run_button(label="::lucide:trash-2:: Delete selected rows")
    mo.vstack([ delete_button, current_view])
    return current_view, delete_button


@app.cell(hide_code=True)
def _(genome, get_tracks, mo, session_template):
    # CSV download using pandas with lazy loading
    async def get_xml_data():
        output = session_template.render(
            genome=f"{genome.value}",
            tracks=get_tracks().to_dict(orient="records")
        )
        return output.encode('utf-8')

    xml_download_lazy = mo.download(
        data=get_xml_data,
        filename="igv_session.xml",
        mimetype="text/xml",
        label="Download session file",
    )
    xml_download_lazy
    return get_xml_data, xml_download_lazy


if __name__ == "__main__":
    app.run()
