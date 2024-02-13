import plotly.express as px


if __name__ == "__main__":
    fig = px.scatter(x=range(10), y=range(10))
    fig.write_html("index.html")
