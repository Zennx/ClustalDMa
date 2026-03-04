/**
 * Custom Molstar Preset: Protein-Nucleic SASA Visualization
 * 
 * Uses Molstar's plugin state API to build visualization step-by-step
 */

async function applyProteinNucleicSASAPreset(plugin, structureRef, params = {}) {
    const {
        showSurface = true,
        surfaceOpacity = 0.85,
        sasaMin = 0,
        sasaMax = 100,
        colorBySASA = true
    } = params;

    console.log('=== Building Molstar state ===');
    console.log('Params:', params);

    // Get the structure ref
    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure) {
        console.error('No structure loaded');
        return;
    }

    const structureCell = structure.cell;
    console.log('Structure cell:', structureCell);

    // Build state update to clear and rebuild
    const update = plugin.build();
    
    // Clear existing components
    const components = [...plugin.managers.structure.hierarchy.current.components];
    console.log(`Clearing ${components.length} components...`);
    for (const comp of components) {
        update.delete(comp.cell);
    }

    // Create protein component (chains A, B) using query expression
    console.log('Creating protein component (chains A, B)...');
    const proteinComp = update.to(structureCell)
        .apply(plugin.builders.structure.component.StructureComponent, {
            type: { 
                name: 'static', 
                params: { } 
            },
            nullIfEmpty: true,
            label: 'Protein'
        }, { 
            tags: 'component',
            ref: 'protein-comp'
        })
        .apply(plugin.builders.structure.component.StructureComponentFromExpression, {
            expression: {
                label: 'Protein Chains',
                expression: plugin.builders.structure.expression.core.rel.eq([
                    plugin.builders.structure.expression.struct.atomProperty.macromolecular.label_asym_id(),
                    plugin.builders.structure.expression.core.type.set(['A', 'B'])
                ])
            }
        }, { ref: 'protein-expr' });

    // Define protein color theme
    let colorTheme;
    if (colorBySASA) {
        console.log(`SASA coloring: ${sasaMin} → ${sasaMax}`);
        colorTheme = {
            name: 'uncertainty',
            params: {
                min: sasaMin,
                max: sasaMax
            }
        };
    } else {
        colorTheme = { name: 'chain-id', params: {} };
    }

    // Add protein representation
    if (showSurface) {
        console.log('Adding molecular surface...');
        proteinComp.apply(plugin.builders.structure.representation.StructureRepresentation3D, {
            type: { name: 'molecular-surface', params: { alpha: surfaceOpacity } },
            colorTheme: colorTheme
        });
    } else {
        console.log('Adding cartoon...');
        proteinComp.apply(plugin.builders.structure.representation.StructureRepresentation3D, {
            type: { name: 'cartoon', params: {} },
            colorTheme: colorTheme
        });
    }

    // Create nucleic acid component (chain L)
    console.log('Creating nucleic acid component (chain L)...');
    const nucleicComp = update.to(structureCell)
        .apply(plugin.builders.structure.component.StructureComponent, {
            type: { 
                name: 'static', 
                params: { } 
            },
            nullIfEmpty: true,
            label: 'Nucleic Acid'
        }, { 
            tags: 'component',
            ref: 'nucleic-comp'
        })
        .apply(plugin.builders.structure.component.StructureComponentFromExpression, {
            expression: {
                label: 'Nucleic Chain',
                expression: plugin.builders.structure.expression.core.rel.eq([
                    plugin.builders.structure.expression.struct.atomProperty.macromolecular.label_asym_id(),
                    plugin.builders.structure.expression.core.type.str('L')
                ])
            }
        }, { ref: 'nucleic-expr' });

    // Add nucleic acid representation
    nucleicComp.apply(plugin.builders.structure.representation.StructureRepresentation3D, {
        type: { name: 'ball-and-stick', params: { sizeFactor: 0.3 } },
        colorTheme: { name: 'element-symbol', params: {} }
    });

    // Commit the state update
    console.log('Committing state update...');
    await update.commit();

    // Reset camera
    plugin.managers.camera.reset();

    console.log('=== State built successfully ===');
}

// Export for use in HTML
if (typeof window !== 'undefined') {
    window.applyProteinNucleicSASAPreset = applyProteinNucleicSASAPreset;
}
