from docutils import nodes

def setup(app):
    app.add_role('filelink', autolink('https://github.com/altMITgcm/MITgcm/blob/master/%s'))
    app.add_role('varlink', autolink('http://mitgcm.org/lxr/ident/MITgcm?_i=%s'))

def autolink(pattern):
    def role(name, rawtext, text, lineno, inliner, options={}, content=[]):
        url = pattern % (text,)
        node = nodes.reference(rawtext, text, refuri=url, **options)
        return [node], []
    return role
