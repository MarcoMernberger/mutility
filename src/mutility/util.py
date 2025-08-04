from typing import Iterable
from scipy.spatial.distance import hamming


def download_file(url, file_object):
    """Download an url"""
    if isinstance(file_object, (str, Path)):
        raise ValueError("download_file needs a file-object not a name")

    try:
        if url.startswith("ftp"):
            return download_ftp(url, file_object)
        else:
            return download_http(url, file_object)
    except Exception as e:
        raise ValueError("Could not download %s, exception: %s" % (repr(url), e))


def download_http(url, file_object):
    """Download a file from http"""
    import requests
    import shutil

    r = requests.get(url, stream=True)
    if r.status_code != 200:
        raise ValueError("HTTP Error return: %i fetching %s" % (r.status_code, url))
    r.raw.decode_content = True
    shutil.copyfileobj(r.raw, file_object)


def download_ftp(url, file_object):
    """Download a file from ftp"""
    import ftplib
    import urllib

    schema, host, path, parameters, query, fragment = urllib.parse.urlparse(url)
    if "ftp_proxy" in os.environ:
        tf = tempfile.NamedTemporaryFile()
        cmd = ["wget", "-O", tf.name, url]
        subprocess.check_call(cmd)
        tf.seek(0, 0)
        bs = 1024 * 1024 * 10
        block = tf.read(bs)
        while block:
            file_object.write(block)
            block = tf.read(bs)
        tf.close()

    else:
        with ftplib.FTP(host) as ftp:
            try:
                ftp.set_pasv(True)
                ftp.login("anonymous", "")
                if "\n" in path:  # pragma: no cover
                    raise ValueError("New line in path: %s" % (repr(path),))
                if path.endswith("/"):
                    ftp.retrbinary("LIST " + path, file_object.write)
                else:
                    ftp.retrbinary("RETR " + path, file_object.write)
            except ftplib.Error as e:
                raise ValueError("Error retrieving urls %s: %s" % (url, e))


def flatten(list_of_lists: Iterable) -> Iterable:
    """
    Flattens a nested list
    """
    if not hasattr(list_of_lists, "__iter__"):
        return [list_of_lists]
    flat = []
    while len(list_of_lists) != 0:
        item = list_of_lists.pop(0)
        if hasattr(item, "__iter__"):
            flat.extend(item)
        else:
            flat.append(item)
    return flat


def hamming(str1, str2):
    return hamming(list(str1), list(str2)) * len(str1)
