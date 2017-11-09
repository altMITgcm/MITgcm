from docutils import nodes, utils

from sphinx.util.nodes import split_explicit_title


def setup(app):
    app.add_role('filelink', filelink('https://github.com/altMITgcm/MITgcm/blob/master/%s'))
    app.add_role('varlink', autolink('http://mitgcm.org/lxr/ident/MITgcm?_i=%s'))

def filelink(pattern):
    def role(name, rawtext, text, lineno, inliner, options={}, content=[]):
        text = utils.unescape(text)
        has_explicit_title, title, part = split_explicit_title(text)
        if not has_explicit_title and part[:1] == '~':
            part = part[1:]
            title = part.split('/')[-1]
        url = pattern % (part,)
        node = nodes.reference(rawtext, title, refuri=url, **options)
        return [node], []
    return role

def autolink(pattern):
    def role(name, rawtext, text, lineno, inliner, options={}, content=[]):
        url = pattern % (text,)
        node = nodes.reference(rawtext, text, refuri=url, **options)
        return [node], []
    return role
