# Import Django views functions & forms.
from django.shortcuts import render, redirect
from .forms import EllinghamForm, OxideForm

# Import GibbsML model.
from .gibbsml.ellingham.ellingham_backend import Ellingham

# Import pymatgen functions.
from pymatgen.core import Composition
from pymatgen.analysis.reaction_calculator import BalancedReaction
from pymatgen.ext.matproj import MPRester
from pymatgen.util.string import unicodeify

# Import Plotly functions.
from plotly.offline import plot
import plotly.graph_objects as go

# Import general functions.
import textwrap
import numpy as np
import math

# Create your views here.

def ellingham(request):

	# Retrieve the dynamic list of metal oxides. 
	all_oxides = request.GET.get('current_oxides', '')
	if all_oxides == '': # If there are not yet any oxides in the list:
		all_oxides = []
	else:
		all_oxides = all_oxides.strip().split(' ')

	# Compile the data (oxide list & blank forms) to be sent to the webpage template.
	context = {
		'all_oxides':all_oxides,
		'oxide_form':OxideForm(),
		'form':EllinghamForm()
	}

	# When the oxide entry form is submitted:
	if request.method == 'GET':

		# String for errors found during oxide input:
		oxide_error = ''

		# Retrieve the oxide field input.
		new_oxide = request.GET.get('id_mo', None)
		if new_oxide is not None:
			new_oxide = new_oxide.strip()

			# Validate & handle errors for oxide entries.
			if new_oxide == '':
				oxide_error = "Please enter a metal oxide before submitting."
			elif ' ' in new_oxide:
				oxide_error = "The metal oxide entry may not contain spaces."
			elif new_oxide in all_oxides:
				oxide_error = "'%s' has already been added to the list!" %(new_oxide)
			elif new_oxide[:3] == 'mp-' and new_oxide.split('mp-', 1)[1].isdigit():
				pass # Allow entries with Materials Project ID syntax to go through.
			else:
				try:
					comp = Composition(new_oxide)
					if comp.is_element:
						oxide_error = "'%s' is an element!" %(new_oxide)
				except:
					oxide_error = "'%s' is an invalid formula!" %(new_oxide)

			if oxide_error != '':
				context.update({'oxide_error':oxide_error})
				# Send the error to the webpage template.
				return render(request, 'ellingham/form.html', context)
			else:
				all_oxides.append(new_oxide)
			
			# Update the dynamic oxide list.
			context['all_oxides'] = all_oxides

	# When the calculation form is submitted:
	if request.method == 'POST':

		# Retrieve the values entered into the form.
		filled_form = EllinghamForm(request.POST or None)
		context['form'] = filled_form

		# Retrieve the dynamic list of metal oxides.
		all_oxides = request.POST.get('current_oxides', '')
		all_oxides = all_oxides.strip().split(' ')
		context['all_oxides'] = all_oxides

		# Retrieve the selected oxides from the dynamic list.
		select_oxides = request.POST.getlist('select_oxides')
		context.update({'select_oxides':select_oxides})

		if filled_form.is_valid(): # If all of the fields have been filled out:

			# Lists to store oxide info:
			plaintext_formulae = []
			presentable_formulae = []
			table_ids = []

			# List to store invalid oxides found during calculation:
			oxide_failures = []

			# Variables to identify duplicate calculations:
			oxide_set = set()
			oxide_duplicates = []

			# String for errors found during calculation:
			calc_error = ''

			# Strings for announcing successful calculations (with unicode & without):
			select_string = ''
			plaintext_string = ''

			# List for the calculated data:
			output = []

			# Retrieve values from the completed form.
			USER_API_KEY = filled_form.cleaned_data.get('USER_API_KEY')
			T = filled_form.cleaned_data.get('T')

			# Catch invalid metal oxides & API keys, and gather info for the valid ones. 
			with MPRester(USER_API_KEY) as mpr:
				for oxide in select_oxides.copy():
					try:
						plaintext_formulae.append(mpr.get_data(oxide, prop='pretty_formula')[0]['pretty_formula'])
					except IndexError:
						oxide_failures.append(oxide)
						select_oxides.remove(oxide)
						all_oxides.remove(oxide)
					except:
						api_key_error = "Invalid Materials Project API Key. Please try again."
						context.update({'api_key_error':api_key_error})
						# Send the error to the webpage template.
						return render(request, 'ellingham/form.html', context)
				for idx, oxide in enumerate(select_oxides):
					if oxide[:3] == 'mp-':
						presentable_formulae.append(unicodeify(plaintext_formulae[idx]))
						table_ids.append(oxide)
					else:
						presentable_formulae.append(unicodeify(oxide))
						table_ids.append(sorted(mpr.get_data(oxide), key=lambda k: k['e_above_hull'])[0]['material_id'])

			# Catch duplicate metal oxides by Materials Project ID.
			for idx, oxide in enumerate(table_ids.copy()):
				if oxide in oxide_set:
					oxide_duplicates.append(select_oxides[idx])
					all_oxides.remove(select_oxides[idx])
					select_oxides[idx] = ''
					plaintext_formulae[idx] = ''
					presentable_formulae[idx] = ''
					table_ids[idx] = ''
				else:
					oxide_set.add(oxide)

			select_oxides[:] = [x for x in select_oxides if x != '']
			plaintext_formulae[:] = [x for x in plaintext_formulae if x != '']
			presentable_formulae[:] = [x for x in presentable_formulae if x != '']
			table_ids[:] = [x for x in table_ids if x != '']

			# Update the dynamic oxide list.
			context['all_oxides'] = all_oxides

			# Set up the base plot.
			def get_dG0_CO_expt(T):
				dG0_CO_expt = -2.23623129e+02 -1.75809004e-01 * T  # Experimental free-energy of CO from the NIST-JANAF database. Malcolm W. Chase, Jr. NIST-JANAF Thermochemical Tables. Washington, DC : New York :American Chemical Society ; American Institute of Physics for the National Institute of Standards and Technology, 1998.
				return dG0_CO_expt

			T_range = np.arange(0, 3001, 1) # Temperature range, from 0K to 3000K.
			dG0_CO = []
			for i in T_range:
				dG0_CO.append(get_dG0_CO_expt(T=i))

			fig = go.Figure()
			fig.add_trace(
				go.Scatter(
					x=T_range,
					y=dG0_CO,
					name="2 C(s) + O\N{SUBSCRIPT TWO}(g) \N{RIGHTWARDS ARROW} 2 CO(g)",
					line=dict(
						color='black',
						width=4,
						dash='dot'
					)
				)
			)

			# Run the calculations.
			for idx, id_mo in enumerate(select_oxides.copy()):
				try:
					ellingham = Ellingham(
						USER_API_KEY=USER_API_KEY,
						id_mo=id_mo
					)

					pred_dH0 = ellingham.get_dH0()
					pred_dS0 = ellingham.get_dS0()
					pred_dG0 = ellingham.get_dG0(T=T)

					# Format the chemical reaction equations.
					try:
						plain_rxn = BalancedReaction.from_string(ellingham.get_balanced_reaction()).normalized_repr_and_factor()[0]
						broken_rxn = plain_rxn.split()
						for i, x in enumerate(broken_rxn):
							if x != '+' and x != '->' and not x.isdigit():
								broken_rxn[i] = unicodeify(x)
						rxn = ' '.join(broken_rxn)
						rxn_sides = rxn.split('->')
						rxn = u'%s \u2192 %s' %(rxn_sides[0], rxn_sides[1])
					
					except:
						plain_rxn = ellingham.get_balanced_reaction()
						broken_rxn = plain_rxn.split()
						for i, x in enumerate(broken_rxn):
							if x != '+' and x != '-->' and not '.' in x:
								broken_rxn[i] = unicodeify(x)
						rxn = ' '.join(broken_rxn)
						rxn_sides = rxn.split('-->')
						rxn = u'%s \u2192 %s' %(rxn_sides[0], rxn_sides[1])

					# Add all calculation data to the output table.
					output.append([presentable_formulae[idx], plaintext_formulae[idx], table_ids[idx], round(pred_dH0, 3), round(pred_dS0, 3), round(pred_dG0, 3), rxn, plain_rxn])

					# Add data for each metal oxide calculation to the output diagram.
					dG0_MO = []
					for i in T_range:
						dG0_MO.append(ellingham.get_dG0(T=i))
					
					fig.add_trace(
						go.Scatter(
							x=T_range,
							y=dG0_MO,
							name="%s [%s]" %(rxn, table_ids[idx]),
							line=dict(
								width=4
							)
						)
					)

				except:
					oxide_failures.append(id_mo)
					select_oxides.remove(id_mo)
					all_oxides.remove(id_mo)
					plaintext_formulae[idx] = ''
					presentable_formulae[idx] = ''
					table_ids[idx] = ''

			plaintext_formulae[:] = [x for x in plaintext_formulae if x != '']
			presentable_formulae[:] = [x for x in presentable_formulae if x != '']
			table_ids[:] = [x for x in table_ids if x != '']

			# Construct the error string for failed oxide calculations (handled for multiple).
			fail_length = len(oxide_failures)
			if fail_length > 2:
				for failure in oxide_failures[:-1]:
					calc_error += "'%s', " %(failure)
				calc_error += "& '%s' are all invalid metal oxides and have been removed from your list. " %(oxide_failures[-1])
			elif fail_length == 2:
				calc_error += "'%s' & '%s' are both invalid metal oxides and have been removed from your list. " %(oxide_failures[0], oxide_failures[1])
			elif fail_length == 1:
				calc_error += "'%s' is an invalid metal oxide and has been removed from your list. " %(oxide_failures[0])

			# Add oxide duplicates to the error string (handled for multiple).
			duplicate_length = len(oxide_duplicates)
			if duplicate_length > 2:
				for duplicate in oxide_duplicates[:-1]:
					calc_error += "'%s', " %(duplicate)
				calc_error += "& '%s' were all duplicates in your metal oxide list, so they have been removed. " %(oxide_duplicates[-1])
			elif duplicate_length == 2:
				calc_error += "'%s' & '%s' were both duplicates in your metal oxide list, so they have been removed. " %(oxide_duplicates[0], oxide_duplicates[1])
			elif duplicate_length == 1:
				calc_error += "'%s' was a duplicate in your metal oxide list, so it's been removed. " %(oxide_duplicates[0])
			
			# Add the error string to the data for the webpage if there were calculation errors.
			if calc_error != '':
				context.update({'calc_error':calc_error})

			# Build the string of calculations (handled for multiple).
			select_length = len(select_oxides)
			if select_length > 2:
				for idx, formula in enumerate(presentable_formulae[:-1]):
					select_string += "%s [%s], " %(formula, table_ids[idx])
					plaintext_string += "%s [%s], " %(plaintext_formulae[idx], table_ids[idx])
				select_string += "& %s [%s]" %(presentable_formulae[-1], table_ids[-1])
				plaintext_string += "& %s [%s]" %(plaintext_formulae[-1], table_ids[-1])
			elif select_length == 2:
				select_string += "%s [%s] & %s [%s]" %(presentable_formulae[0], table_ids[0], presentable_formulae[1], table_ids[1])
				plaintext_string += "%s [%s] & %s [%s]" %(plaintext_formulae[0], table_ids[0], plaintext_formulae[1], table_ids[1])
			elif select_length == 1:
				select_string += "%s [%s]" %(presentable_formulae[0], table_ids[0])
				plaintext_string += "%s [%s]" %(plaintext_formulae[0], table_ids[0])
			else: # If none of the calculations succeed:
				calc_error += "None of the metal oxides you selected were valid; please try again."
				context['calc_error'] = calc_error
				# Send the error to the webpage template.
				return render(request, 'ellingham/form.html', context)

			# Make the titles for the output diagram, table, & downloadable table.
			message = "Results for %s (%s K)" %(select_string, T)
			plaintext_message = "Results for %s (%s K)" %(plaintext_string, T)

			# Make the headers for the output table.
			oxide_header = "Metal Oxide"
			id_header = "Materials Project ID"
			formation_header = "Formation Energy (kJ/mol)"
			entropy_header = "dG/dT (kJ/(mol*K))"
			free_header = "Free Energy (kJ/mol)"
			rxn_header = "Balanced Reaction Equation"
			headers = [oxide_header, id_header, formation_header, entropy_header, free_header, rxn_header]

			# Format the output diagram.
			split_title = textwrap.wrap(select_string, width=40)
			fig.update_layout(
				height=750,
				template='plotly_white',
				font_family='sans-serif',
				font=dict(
					size=16
				),
				title="Ellingham Diagram for %s"  %('<br>'.join(split_title)),
				xaxis_title="Temperature (K)",
				yaxis_title="Free Energy (kJ/mol)",
				legend=dict(
					traceorder='reversed',
					font_size=14
				)
			)

			# Prepare diagram for HTML rendering.
			diagram = plot(
				fig,
				output_type='div',
				include_plotlyjs=False
			)

			# Add new attributes to the compiled webpage template data.
			context.update({
				'message':message,
				'plaintext_message':plaintext_message,
				'headers':headers,
				'output':output,
				'diagram':diagram
			})
	
	# Send all relevant information to the webpage template & render it.
	return render(request, 'ellingham/form.html', context)