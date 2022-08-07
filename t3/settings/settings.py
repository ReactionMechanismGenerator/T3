"""
T3's settings

You may keep a short version of this file in a local ".t3" folder under your home folder.
Any definitions made to the local file will take precedence over this file.
"""


# The execution type can be either 'incore', i.e., executed in the same processor,
# or 'local', i.e., to be submitted to the server queue if running on a server.
# If running on a local server, ARC's settings for ``local`` will be used.
execution_type = {
    'rmg': 'incore',
    'arc': 'incore',
}
