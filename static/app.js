const form = document.getElementById('analysis-form');
const result = document.getElementById('result');
const errorBox = document.getElementById('error');

form.addEventListener('submit', async (event) => {
  event.preventDefault();
  result.classList.add('hidden');
  errorBox.classList.add('hidden');

  const payload = {
    sequence: document.getElementById('sequence').value,
    frame: Number(document.getElementById('frame').value),
    reverse_complement: document.getElementById('reverse_complement').checked,
    stop_at_stop: document.getElementById('stop_at_stop').checked,
    motif: document.getElementById('motif').value || null,
  };

  try {
    const response = await fetch('/api/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });

    const data = await response.json();
    if (!response.ok) {
      throw new Error(data.detail || 'Analysis failed');
    }

    const motifText = data.motif
      ? data.motif_positions.length
        ? `Motif '${data.motif}' positions (0-based): ${data.motif_positions.join(', ')}`
        : `Motif '${data.motif}' not found`
      : 'No motif requested';

    result.innerHTML = `
      <h2>Results</h2>
      <p><strong>DNA:</strong> ${data.dna}</p>
      <p><strong>Protein:</strong> ${data.protein}</p>
      <p><strong>Frame:</strong> ${data.frame}</p>
      <p><strong>Reverse complement:</strong> ${data.reverse_complement ? 'yes' : 'no'}</p>
      <p><strong>Stop at stop codon:</strong> ${data.stop_at_stop ? 'yes' : 'no'}</p>
      <p><strong>Motif:</strong> ${motifText}</p>
    `;
    result.classList.remove('hidden');
  } catch (error) {
    errorBox.textContent = error.message;
    errorBox.classList.remove('hidden');
  }
});
