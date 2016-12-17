#!/usr/bin/env python
'''
Created on 13 Sep 2011

@author: Fergal
'''
import cgitb; cgitb.enable()
import cgi
import sys
import tempfile

#My modules
import scramble_pep

page_template = '''
<html>
    <head>
        <title>Peptide Controls</title>
        <script type="text/javascript" src="/~testing/biowareweb/Javascript/sorttable.js"></script>
        <script type="text/javascript" src="../control_pep_website/control_peptides.js"></script>
        <link rel="stylesheet" type="text/css" href="/~testing/biowareweb/Stylesheets/bioware.css">
        <link rel="stylesheet" type="text/css" href="/~testing/biowareweb/Stylesheets/example.css">
        <style type="text/css">.switchcontent{display:none;}
            table.sortable a.sortheader {
            background-color:#eee;
            color:#666666;
            font-weight: bold;
            text-decoration: none;
            display: block;
            }
            table.sortable span.sortarrow {
            color: black;
            text-decoration: none;
            }


            table.sortable td 
            {
            white-space: nowrap;
            font-size:7pt;
            background-color:#F8F8F8;
            bordercolor:black; 
            }
            
            table.sortable th 
            {
            white-space: nowrap;
            font-size:9pt;
            }
        </style>
        
    </head>
    <body >
        <h1>PepControls </h1>
        <br />
        <br />
        <h2>Control Peptide Output <div id="sequence" style = "display:none">%(sequence)s</div></h2>
        <div>Sort peptide controls by clicking on column headers. <br /><br />
        <a href='javascript:assemble_download_query()'>Download control peptides as a .csv file. </a><br /><br />
        <a href='../control_pep_website/worked_example.html'>Choosing Control Peptides: Worked example</a><br /><br /><br />
        </div>
            %(table)s
    <body>
</html>
'''

def parse_form():
    form = cgi.FieldStorage()
    sequence = form.getvalue('sequence')
    download = form.getvalue('download')
    return sequence, download
    
    
    
def html_tabulate(data):
    headers_written = False
    table = '<table width=100% id="main_table" class="sortable">\n'
    for row in data:
        table += "<tr>\n"
        if not headers_written:
            for item in row:
                table += "<th>" + str(item) + "</th>"
            headers_written = True
        else:
            for item in row:
                if item != 0:
                    item = ' ' if not item else item
                table += "<td>" + str(item) + "</td>"
        table += "</tr>"
    table += "</table>"
    return table

def main():
    html_string = "Content-type: text/html\n\n"
    error_string = "Content-type: text/plain\n\n Error!"
    download_string = "Content-type: text/csv\nContent-disposition:attachment;filename=%s.csv\n\n"
    sequence, download = parse_form()
    # try:
    if download:
        data = scramble_pep.write_controls(sequence, outfile=False, chatty=False)
        sys.stdout.write(download_string % sequence)
        sys.stdout.write(data)
    else:
        rows = scramble_pep.yield_controls(sequence, chatty=False)
        table = html_tabulate(rows)
        sys.stdout.write(html_string)
        sys.stdout.write(page_template % {'table':table, 'sequence':sequence})
    # except Exception:
        # sys.stdout.write(error_string)
        # sys.stdout.write(sequence)

if __name__ == '__main__':
    main()
