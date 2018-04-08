import numpy as np
import time
from vispy import app
from vispy import gloo

c = app.Canvas(keys='interactive')

vertex = """
attribute vec2 a_position;
void main (void)
{
    gl_Position = vec4(a_position, 0.0, 1.0);
}
"""
fragment = """
void main()
{
    gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
}
"""

program = gloo.Program(vertex, fragment)
program['a_position'] = np.c_[
    [-0.3, 0.3],
    [0.5, -0.5]
].astype(np.float32)

@c.connect
def on_resize(event):
    gloo.set_viewport(0, 0, *event.size)

@c.connect
def on_draw(event):
    gloo.clear((1, 1, 1, 1))
    program.draw('line_Strip')

c.show()
app.run()

for i in range(0, 10):
    time.sleep(0.1)
    program['a_position'] = np.c_[
                [-0.3, 0.3],
                [0.5 - i / 10.0, -0.5 + i / 10.0]
            
            ].astype(np.float32)
    gloo.clear((1, 1, 1, 1))
    program.draw('line_Strip')
