from shiny import reactive
from shiny.express import input, render, ui

ui.input_action_link("action_link", "Run")  

a = reactive.value(0)

@reactive.effect
@reactive.event(input.action_link)
def set_default(default_value):
    a.set(default_value)

@render.text()
def counter():
    return f"{a()}"