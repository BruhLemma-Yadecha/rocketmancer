const PALETTE = ['red', 'orange', 'yellow', 'green', 'blue', 'purple'];

export function color(index) {
  return PALETTE[index % PALETTE.length];
}

function separateThousands(s) {
  const [int, dec] = s.split('.');
  const grouped = int.replace(/\B(?=(\d{3})+(?!\d))/g, ',');
  return dec !== undefined ? `${grouped}.${dec}` : grouped;
}

export function fmt(n, decimals = 2) {
  if (typeof n !== 'number') return '-';
  return separateThousands(n.toFixed(decimals));
}
