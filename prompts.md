this is a basic unsorted record of interesting lessons/examples from prompts

tips:

1. prefix prompts with `[check llm_ref.md for context before answering questions about this repo; warn if not accessibly and ask if continue]`
    - e.g. `[above prompt prefix]` `[rest of prompt]`
    - using github copilot, the above prefix achieves an incredibly better answer to a question like `'how can volco use spheres to generate a continuous line of material?'`


# Example prompt to generate a llm reference document for a software repo: 

```
Create a concise LLM Reference Document (llm_ref.md) for this repository that will help large language models understand the codebase structure without needing to upload all files as context. The document should be optimized to minimize context needed for LLM code editing while providing sufficient information to understand the system architecture.

Follow this structure:

## 1. Repository Overview
Provide a brief description of the repository's purpose and main functionality (2-3 sentences).

Include a brief introduction that:
- States which files were examined to create the document
- Notes that the document should be updated as the repository evolves
- Mentions any limitations of the document's coverage
- Briefly acknowledges any coding conventions or patterns that are imperative to understanding the repository structure

## 2. Directory Structure
Create a tree representation of the main directories and their purposes. Focus on the most important directories and files.

## 3. Core Components
For each major component/module:
- Use a clear subsection header that includes the component name AND relevant file paths
- Describe its purpose in 1-2 sentences
- List key classes/functions with very brief descriptions
- Do not include exhaustive method descriptions - focus on the most important ones
- Organize components by functional area for better clarity

## 4. Key Data Flows
Describe 2-4 main data flows through the system:
- Use arrows (→) to show how data moves between components
- Include specific method calls where relevant
- Show the complete path from input to output with numbered steps
- Highlight any conditional branches in the flow

## 5. Configuration Options
List configuration parameters with brief descriptions of each. Group related parameters and indicate their purpose.

## 6. Usage Patterns
Show basic usage examples in the most concise format possible (command line, API calls, etc.). Include complete but minimal examples.

## 7. Key Algorithms
Use bullet points to describe the main algorithms in 1-2 lines each, focusing on their purpose and approach.

## 8. Module Dependencies
Show the main dependencies between modules using arrows (→). Focus on high-level relationships.

## 9. Common Modification Patterns
List 2-4 common ways the code might be modified or extended, with clear step-by-step instructions for each.

## 10. Performance Considerations
List key factors affecting performance in bullet points with brief explanations of their impact.

Important guidelines:
- Begin with a clear introduction paragraph explaining the purpose of the document
- Include a prominent note about which files were examined and any limitations of coverage
- Be extremely concise - prefer bullet points over paragraphs - don't explain functions with self-explanatory titles
- Focus on relationships between components rather than implementation details
- Include file paths with component names to help locate relevant code
- Prioritize information that helps understand the system architecture
- Try to capture the entire hierarchical flow of function calls for major functionality
- Avoid redundancy - information should appear in only one section
- Use hierarchical organization with clear subsection headers for better navigation
- The final document should be readable in under 5 minutes

Examine key files that represent the core functionality, focusing on main entry points, core modules, and configuration files. Acknowledge which files were checked.
```