from django import forms

__author__ = "Tobias Eegholm Hoffmann, Alexander Urban"
__email__ = "aurban@atomistic.net"
__date__ = "2021-03-18"
__version__ = "1.0"


class OxideForm(forms.Form):
    id_mo = forms.CharField(
        label="Metal Oxide Input",
        max_length=50,
        help_text="chemical formula or Materials Project ID",
        widget=forms.TextInput(
            attrs={
                'class': "form-control",
                'aria-describedby': "id_id_mo_help_text",
                'placeholder': "e.g. LiAlO2 or mp-3427"
            }
        )
    )


class EllinghamForm(forms.Form):
    USER_API_KEY = forms.CharField(
        label="Materials Project API Key",
        max_length=20,
        help_text="Don’t have one? Register for free at <a href=\"https://materialsproject.org\" target=\"_blank\">materialsproject.org</a>.",
        widget=forms.TextInput(
            attrs={
                'class': "form-control",
                'aria-describedby': "id_USER_API_KEY_help_text"
            }
        )
    )
    T = forms.IntegerField(
        label="Temperature (K)",
        max_value=3000,
        min_value=0,
        initial=273,
        help_text="for the predicted free energy calculation",
        widget=forms.NumberInput(
            attrs={
                'class': "form-control",
                'aria-describedby': "id_T_help_text"
            }
        )
    )
