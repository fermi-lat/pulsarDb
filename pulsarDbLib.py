#$Id$
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['pulsarDb'])
    env.Tool('st_appLib')
    env.Tool('st_facilitiesLib')
    env.Tool('st_streamLib')
    env.Tool('timeSystemLib')
    env.Tool('tipLib')

def exists(env):
    return 1
