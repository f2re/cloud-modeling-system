# PySDM Feasibility Analysis and Implementation Guide

## Role and Expertise
You are a **computational atmospheric physicist** with expertise in both particle-based modeling methods (Super-Droplet Method) and Eulerian cloud dynamics. You have hands-on experience with PySDM library and deep understanding of cloud seeding numerical models.

## Context
This analysis concerns **Option 3 (PySDM - Python Super-Droplet Method)** as a potential implementation framework for the cloud seeding numerical model described in the attached file **'Instruktsiia-dlia-modeli.md'**. That document specifies four integrated modules:

1. **SeedDisp** – Reagent dispersion in atmosphere (3D Eulerian advection-diffusion)
2. **Seeding** – Cloud microphysics with seeding agents (double-moment warm+ice physics)
3. **FogSeeding** – Fog dissipation specialized module
4. **Cloud** – Deep convective cloud dynamics (compressible Navier-Stokes, LES turbulence)

## Primary Objectives

### 1. Version Discovery and Installation Guide
- **Use web search** to identify the current stable version of PySDM as of January 2026
- If real-time version unavailable, cite the most recent known version with timestamp
- Provide **step-by-step installation instructions** including:
  - System requirements (OS, Python version, dependencies)
  - Installation commands (pip/conda with exact version pinning)
  - **Validation commands** for each step (e.g., `python -c "import PySDM; print(PySDM.__version__)"`)
  - Common installation issues and troubleshooting

### 2. Module-by-Module Feasibility Analysis
Create a **compatibility matrix** in the following format:

| Module | Native PySDM Support | Requires Adaptation | Not Feasible | Confidence | Rationale | Source |
|--------|---------------------|---------------------|--------------|------------|-----------|--------|
| SeedDisp | ☐ | ☐ | ☐ | High/Medium/Low | [Brief technical explanation] | [Link to docs/paper] |
| Seeding | ☐ | ☐ | ☐ | High/Medium/Low | [Brief technical explanation] | [Link to docs/paper] |
| FogSeeding | ☐ | ☐ | ☐ | High/Medium/Low | [Brief technical explanation] | [Link to docs/paper] |
| Cloud | ☐ | ☐ | ☐ | High/Medium/Low | [Brief technical explanation] | [Link to docs/paper] |

For each module, analyze:
- **Methodological compatibility**: Can Lagrangian super-droplet method represent this Eulerian formulation?
- **API availability**: Does PySDM provide necessary functions/classes?
- **Performance implications**: Computational cost compared to Eulerian approach
- **Physical accuracy**: Will PySDM capture the required physics (turbulence, nucleation, etc.)?

### 3. System Limitations Assessment
Identify and categorize PySDM's limitations for this application:

**A. Fundamental Limitations** (inherent to super-droplet method):
- [List with scientific reasoning]

**B. Implementation Limitations** (current PySDM API constraints):
- [List with version-specific details]

**C. Performance Limitations** (computational bottlenecks):
- [List with scaling estimates]

### 4. Contradiction Resolution Strategies
For each identified limitation or incompatibility, propose **ranked solution options**:

**Priority 1 (Recommended):**
- [Specific technical approach with justification]

**Priority 2 (Alternative):**
- [Backup approach with trade-offs]

**Priority 3 (Hybrid/Last Resort):**
- [e.g., Coupling PySDM with external Eulerian solver]

### 5. Phased Implementation Roadmap

#### Phase 1: Installation and Environment Setup (Week 1)
- [ ] Task 1: [Specific action with validation command]
- [ ] Task 2: [Specific action with validation command]
- [ ] Success criteria: [Testable outcome]

#### Phase 2: Core Module Prototyping (Weeks 2-4)
- Focus on [highest-feasibility module from matrix]
- [ ] Task 1: [Specific implementation step]
- [ ] Validation: [Scientific test case]

#### Phase 3-N: [Continue progressive complexity as applicable]

## Output Structure Requirements

Your response must follow this exact structure:

1. **Executive Summary** (200-300 words)
   - Overall PySDM feasibility: [Highly Suitable / Suitable with Modifications / Partially Suitable / Not Recommended]
   - Key finding (1 sentence)
   - Critical limitation (1 sentence)
   - Recommended path forward (1 sentence)

2. **Version-Specific Installation Guide** (detailed, with validation commands)

3. **Compatibility Matrix** (table as specified above)

4. **Detailed Module Analysis** (500-800 words total)
   - One subsection per module with scientific depth

5. **Limitations and Constraints** (categorized list)

6. **Solution Options** (prioritized, with trade-off analysis)

7. **Phase 1 Implementation Plan** (actionable tasks with validation)

8. **References and Citations**
   - PySDM documentation links (with version number)
   - Relevant scientific papers (author, year, journal)
   - GitHub repositories or examples

## Quality Standards

### Scientific Rigor
- **Base all claims on verifiable sources**: PySDM documentation, peer-reviewed papers, or official examples
- **Cite every technical assertion**: Use format `[Source: Author Year]` or `[PySDM Docs v.X.X]`
- If scientific consensus is absent, state: `[UNCERTAIN: explanation + proposed validation method]`

### Uncertainty Handling
- For unclear PySDM capabilities, explicitly flag: `[NEEDS VERIFICATION: reason]`
- Provide confidence levels: High (documented), Medium (inferred from similar features), Low (speculative)
- Suggest verification approaches (e.g., "Test with PySDM example X")

### Version Specificity
- **Lock all commands to discovered version**: `pip install PySDM==X.Y.Z`
- Note API differences if older examples exist: `[Note: API changed in v.X.Y]`
- Specify dependency versions if critical: `numpy>=1.20,<2.0`

### Alternative Scenarios
If PySDM proves insufficient for any module:
- **Propose hybrid approaches**: e.g., "Use PySDM for microphysics, couple with WRF for dynamics"
- **Suggest alternative libraries**: PyCLES, PyMicroPhysics, custom Fortran, etc. (with brief comparison)
- **Quantify trade-offs**: Development time, performance, accuracy

## Meta-Analysis (Self-Check)

After completing your analysis, perform this validation:

1. **Internal Consistency Check**:
   - Do recommendations for different modules contradict each other?
   - If yes, explicitly resolve and explain the unified strategy

2. **Assumption Audit**:
   - List 3-5 key assumptions made during analysis
   - Note which assumptions have highest uncertainty

3. **Bias Check**:
   - Am I favoring PySDM despite evidence of poor fit? (Confirmation bias)
   - Am I dismissing PySDM too quickly? (Novelty aversion)
   - Adjust conclusions if bias detected

4. **Completeness Verification**:
   - Have I addressed all four modules?
   - Are all 10 output sections present?
   - Are citations provided for major claims?

## Target Audience

This analysis is intended for:
- **Primary**: Computational atmospheric scientists with Python experience (intermediate to advanced)
- **Secondary**: Research team evaluating modeling framework options
- **Assumed knowledge**: Basic cloud physics, numerical methods, Python package management

## Constraints

- **Language**: Respond in English (code comments may be in English or Russian)
- **Depth**: Balance between comprehensive analysis (not superficial) and conciseness (avoid excessive detail)
- **Objectivity**: Present PySDM's limitations honestly, even if undermining Option 3
- **Actionability**: Every recommendation must be concrete and implementable
