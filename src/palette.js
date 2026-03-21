const PALETTE = ['red', 'orange', 'yellow', 'green', 'blue', 'purple'];

export function color(index) {
  return PALETTE[index % PALETTE.length];
}

export function fmt(n) {
  return typeof n === 'number' ? n.toFixed(4) : '-';
}
