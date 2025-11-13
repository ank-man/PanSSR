# panssrator/report_generator.py
import json
from panssrator import utils

def generate_html_report(marker_stats: dict, output_file: str):
    """
    Generate an HTML report from marker statistics.
    marker_stats: a dictionary containing statistical summaries.
    The report includes interactive charts using echarts (to be embedded in the HTML).
    """
    # For illustration, we produce a simple HTML template.
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>PanSSRAtor Marker Report</title>
        <script src="https://cdn.jsdelivr.net/npm/echarts/dist/echarts.min.js"></script>
        <style>
            body {{ font-family: Arial, sans-serif; }}
            .chart {{ width: 600px; height: 400px; margin: 20px auto; }}
        </style>
    </head>
    <body>
        <h1>PanSSRAtor Marker Report</h1>
        <div id="chart1" class="chart"></div>
        <script>
            var chartDom = document.getElementById('chart1');
            var myChart = echarts.init(chartDom);
            var option;
            option = {{
                title: {{
                    text: 'Marker Distribution'
                }},
                tooltip: {{}},
                legend: {{
                    data: ['Count']
                }},
                xAxis: {{
                    data: {json.dumps(marker_stats.get('x_axis', []))}
                }},
                yAxis: {{}},
                series: [{{
                    name: 'Count',
                    type: 'bar',
                    data: {json.dumps(marker_stats.get('counts', []))}
                }}]
            }};
            myChart.setOption(option);
        </script>
    </body>
    </html>
    """
    with open(output_file, "w") as f:
        f.write(html_template)
    utils.logger.info("HTML report generated at %s", output_file)

if __name__ == '__main__':
    # Test report generation
    stats = {
        "x_axis": ["chr1", "chr2", "chr3"],
        "counts": [120, 80, 150]
    }
    generate_html_report(stats, "report.html")

