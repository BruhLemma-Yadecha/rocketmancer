const GITHUB_URL = 'https://github.com/BruhLemma-Yadecha/rocketmancer-web';

export default function Footer() {
  return (
    <footer>
      <span>&copy; 2026 B.L. Yadecha</span>
      {' · '}
      <a href={GITHUB_URL} target="_blank" rel="noopener noreferrer">
        GitHub
      </a>
    </footer>
  );
}
