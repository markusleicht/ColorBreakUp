# ColorBreakUp

The CBU_MODEL predicts a viewer's perception of Color Break-Up (CBU) on an objectified basis for predefined scenarios. The model scenario always includes (1) a field-sequential color display system presenting digital content and (2) a viewer observing the content. Under certain conditions, such a scenario potentially leads to the occurrence of CBU effects. The CBU effects are computed by the CBU_MODEL in a multi-stage process.

The CBU_LOOP is a supplement to the CBU_MODEL. The loop allows to run the model multiple times. During the runs, the model output is calculated under variation of one input parameter. All other parameters are constant in order to determine the impact of the variable parameter on CBU perception.

Both MATLAB functions were written during the author's doctorate. The associated doctoral thesis refers to model and loop and provides further information on CBU: Leicht, M. (2022). Perception of Color Break-Up [Doctoral thesis]. Ilmenau University of Technology.
