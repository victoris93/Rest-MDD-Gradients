# Discussion
The aim of the present study was to compare macroscale cortical gradients of patients diagnosed with MDD and healthy controls. Along the first gradient, we hypothesized a shift of scores towards the transmodal end of the spectrum in drug-naïve patients. We expected a reverse trend in patients on medication, namely a shift of scores towards the unimodal end of the gradient resulting in more dispersed transmodal regions. We extracted two gradients from readily preprocessed fMRI timeseries provided by the REST-meta-MDD Consortium {cite:p}`yan_reduced_2019` and compared them between healthy controls and depressed patients. 

We demonstrated that depressed drug-naïve subjects had a shift of gradient scores towards the transmodal end of the first gradient in the visual network: 
```{figure} _build/html/_images/Results_37_2.png
:name: res-1-ref
```
<style>
.marginauto {
    margin: 10px auto 20px;
    display: block;
}
</style>
<img class="marginauto" src="figures/yeo_result_1_grad_DN_NC.png" width = 450 alt="centered image" />

This result confirmed our prediction as to expanded transmodal areas: the scores of the regions overlapping with the visual network were shifted towards the transmodal end of the spectrum in drug-naïve patients.  Additionally, this score shift implies higher functional connectivity between the visual network and associative regions. This observation can putatively constitute a neural substrate of the most common symptom of depression – prolonged periods of self-focused rumination and the failure to objectively assess other people’s intentions – since hierarchical abnormalities within sensory regions can result in further distortions occurring in upstream regions all the way to higher-order associative areas. The involvement of the visual network in MDD is consistent with the results reported in {cite:t}`yan_reduced_2019` thereby lending support to the gradient framework as a valid alternative to classical functional connectivity analyses. We therefore suggest further research into the role of the visual network in psychopathology and notably in MDD.

The exploratory analysis revealed a shift in scores of two regions largely overlapping with the DMN towards the extremity occupied by the visual network:

```{figure} _build/html/_images/Results_39_2.png
:name: res-2-ref
```
<img class="marginauto" src="figures/yeo_result_2_grad_DN_NC.png" width = 700 alt="centered image" />

 Bearing in mind that the DMN was attributed transitional values close to zero, one can conclude, that in drug-naïve patients, functional connectivity is higher between right angular gyrus and paracingulate cortex and the rest of the DMN along the second gradient. These observations are in line with the previous research on MDD reporting greater activity and connectivity within the DMN. The third region manifesting a significant difference in the gradient score – left middle temporal gyrus, temporooccipital part – was shifted in the opposite direction in drug-naïve depressed subjects, thus implying greater functional connectivity with the regions attributed negative values and overlapping with multiple resting state networks: the DMN, ventral attention and frontoparietal. The latter two networks occupy the opposite end of the second gradient relative to the visual network which might account for the finding. The fact that multiple networks coalesce within the left temporal middle cortex could account for the shift towards the ventral attention and frontoparietal networks rather than the DMN. 

With regards to subjects on medication, contrary to our prediction, no effect was found in the first gradient. Nevertheless, one cannot exclude the possibility of antidepressants exerting an effect on macroscale cortical hierarchy. That is to say, the absence of a significant effect is a consequence of insufficient sensibility either attributable to the methodology employed here or intrinsic to the dataset we had at our disposal. It is likely that the former is more relevant since, as mentioned above, the power of the sample was insufficient to detect a small effect (Cohen’s d ≤ 0.2).

Nonetheless, a score shift was detected along the second gradient in the right frontal operculum cortex which is part of the ventral attention network:

```{figure} _build/html/_images/Results_43_2.png
:name: res-3-ref
```
<img class="marginauto" src="figures/yeo_result_2_grad_MED_NC.png" width = 250 alt="centered image" />

 This observation points at a greater temporal decoupling of activity within the ventral attention network. Indeed, it was previously demonstrated that the administration of antidepressants reduces overall functional connectivity, including within the ventral attention network. {cite:p}`li_span_2021` However, to establish associations between abnormal connectivity patterns in specific networks and the behavioral manifestation of MDD symptomatology, further research into resting state functional connectivity is needed.

These results support the validity of the gradient-based approach for characterizing large-scale differences in functional organization associated with MDD. However, a number of limitations are worth mentioning. These can be divided into caveats intrinsic to the data and the caveats ensuing from the methodology. First, Harvard-Oxford parcellation is too coarse to be implicated in analyses of resting state functional connectivity, including gradient-based methods. Also, since Yeo’s atlas was not included in the preprocessed dataset, we had to adopt the strategy of visual assessment of spatial overlap between ROI and Yeo’s networks. Second, the dataset used here comprised only Chinese subjects which makes the extension of the results onto other ethnicities a complicated matter (e.g. {cite:t}`tang_brain_2018`). We suggest further investigations based on the gradient approach be conducted on non-Chinese datasets such as UK Biobank and ENIGMA MDD.

As for GCCA, akin to any dimensionality reduction algorithm, it is based on the assumption that the components resulting from the transformation of the data are real entities and are treated as such by researchers. However, it is important to keep in mind that every component is only the product of the captured maximum variance in the dataset which constitutes a minor fraction of the original variance.  In other words, every dimensionality reduction algorithm, be it GCCA, diffusion embedding or some other, excludes most of the original variance from subsequent analyses. Furthermore, the peculiarity of PCA-like algorithms is what can put into question its ecological validity: the resulting components are orthogonal, i.e. uncorrelated. Methodological advancements leading to elimination of these limitations are called for. Also, the question still remains as to the appropriate number of components to be extracted in the case of fMRI data.  Moderate caution in the interpretation of outcomes of dimensionality reduction algorithms is paramount. 

The matter of effect size is also of great importance in the context of neuroimaging studies. The effect we detected in every case was small (Cohen’s d < 0.3). When confronted with a statistically significant but small effect, one may consider the following questions: what is the meaning of a small, but significant difference between two large, independent samples? How informative is this difference? These questions constitute the core of the statistical debate within the field whereby the proponents of reporting effect sizes claim that statistical significance (i.e. p-value) alone does not suffice to deem the result meaningful. {cite:p}`Chen064212` To tackle this issue, more meta-analyses and comparative methodological studies need to be conducted. In this study, we stood by the principle of reporting effect sizes for the sake of facilitation of subsequent meta-analyses and replications.

Our findings pave the way for further investigations based on the gradient approach since it has the potential of offering a holistic overview of alterations in macroscale functional cortical hierarchy in psychopathology. Ultimately, we posit that macroscale cortical gradients may be used to obtain much coveted answers to most fundamental questions in cognitive neuroscience: how are our brains different? How can this difference be explained in terms of developmental mechanisms and, more specifically, neurogenesis? How can we globally characterize biomarkers of psychopathology? Can one elaborate potential scenarios of emergence of psychopathological conditions based on neuroimaging data? These questions could benefit from holistic frameworks such as macroscale cortical gradients.

# REFERENCES

```{bibliography}
```